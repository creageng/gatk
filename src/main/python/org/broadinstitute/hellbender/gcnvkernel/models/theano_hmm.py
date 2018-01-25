import numpy as np
import theano as th
import theano.tensor as tt
import pymc3 as pm
from typing import Optional, Set
from .. import types
from . import commons
from scipy.misc import logsumexp


class TheanoForwardBackward:
    """Implementation of the forward-backward algorithm using `theano.scan`."""
    def __init__(self,
                 log_posterior_output: Optional[types.TensorSharedVariable] = None,
                 admixing_rate: float = 1.0,
                 include_alpha_beta_output: bool = False,
                 resolve_nans: bool = False):
        """Initializer.

        Args:
            log_posterior_output: if not None, the new log posterior will be written to this shared tensor;
                otherwise, it will be returned by `TheanoForwardBackward.perform_forward_backward`
            admixing_rate: a float in range [0, 1] denoting the amount of the new posterior to admix with the
                old posterior (higher = more of the new posterior)
            include_alpha_beta_output: include forward and backward tables in the return values
                of `TheanoForwardBackward.perform_forward_backward`
            resolve_nans: if True, expression such as inf - inf resulting in NaNs will be properly handled
        """
        self.admixing_rate = admixing_rate
        self.include_alpha_beta_output = include_alpha_beta_output
        self.resolve_nans = resolve_nans
        assert 0.0 < admixing_rate <= 1.0, "Admixing rate must be in range (0, 1]"
        self.log_posterior_output = log_posterior_output
        self._forward_backward_theano_func = self._get_compiled_forward_backward_theano_func()

    def perform_forward_backward(self,
                                 num_states: int,
                                 temperature: float,
                                 log_prior_c: np.ndarray,
                                 log_trans_tcc: np.ndarray,
                                 log_emission_tc: np.ndarray,
                                 prev_log_posterior_tc: np.ndarray):
        """Runs the forward-backward algorithm.

        Args:
            num_states: number of states for each node of the Markov chain
            temperature: a scale factor of the chain entropy
            log_prior_c: prior probability vector for the first node
            log_trans_tcc: transition probability matrices for each directed vertex
            log_emission_tc: emission probability vector for each node
            prev_log_posterior_tc: previous estimate of the log posterior (used for admixing)

        Returns:
            see `TheanoForwardBackward._get_compiled_forward_backward_theano_func`
        """
        return self._forward_backward_theano_func(
            num_states, temperature, log_prior_c, log_trans_tcc, log_emission_tc, prev_log_posterior_tc)

    @th.configparser.change_flags(compute_test_value="ignore")
    def _get_compiled_forward_backward_theano_func(self):
        """Returns a compiled theano function that perform forward-backward and either updates log posterior
        probabilities or returns it.

        Note:
            The returned theano function takes 6 inputs:

                num_states (integer scalar),
                temperature (float scalar),
                log_prior_c (float vector),
                og_trans_tcc (float tensor3),
                log_emission_tc (float matrix)
                prev_log_posterior_tc (float matrix)

            If a `log_posterior_output` shared tensor is given to the class initializer,
            the return tuple will be:

                update_norm_t, log_data_likelihood,
                (+ alpha_tc, beta_tc if self.include_alpha_beta_output == True)

            and the posterior will be directly written to `self.log_posterior_output`. Otherwise,
            return tuple will be:

                admixed_log_posterior_tc, update_norm_t, log_data_likelihood,
                (+ alpha_tc, beta_tc if self.include_alpha_beta_output == True)

        Returns:
            A compiled theano function
        """
        num_states = tt.iscalar('num_states')
        temperature = tt.scalar('temperature')
        log_prior_c = tt.vector('log_prior_c')
        log_trans_tcc = tt.tensor3('log_trans_tcc')
        log_emission_tc = tt.matrix('log_emission_tc')
        prev_log_posterior_tc = tt.matrix('prev_log_posterior_tc')

        new_log_posterior_tc, log_data_likelihood_t, alpha_tc, beta_tc = self._get_symbolic_log_posterior(
            num_states, temperature, log_prior_c, log_trans_tcc, log_emission_tc, self.resolve_nans)

        admixed_log_posterior_tc = commons.safe_logaddexp(
            new_log_posterior_tc + np.log(self.admixing_rate),
            prev_log_posterior_tc + np.log(1.0 - self.admixing_rate))

        log_data_likelihood = log_data_likelihood_t[-1]  # in theory, they are all the same
        update_norm_t = commons.get_jensen_shannon_divergence(admixed_log_posterior_tc, prev_log_posterior_tc)

        ext_output = [alpha_tc, beta_tc] if self.include_alpha_beta_output else []
        inputs = [num_states, temperature, log_prior_c, log_trans_tcc, log_emission_tc, prev_log_posterior_tc]
        if self.log_posterior_output is not None:
            return th.function(inputs=inputs,
                               outputs=[update_norm_t, log_data_likelihood] + ext_output,
                               updates=[(self.log_posterior_output, admixed_log_posterior_tc)])
        else:
            return th.function(inputs=inputs,
                               outputs=[admixed_log_posterior_tc, update_norm_t, log_data_likelihood] + ext_output)

    @staticmethod
    def _get_symbolic_log_posterior(num_states: tt.iscalar,
                                    temperature: tt.scalar,
                                    log_prior_c: types.TheanoVector,
                                    log_trans_tcc: types.TheanoTensor3,
                                    log_emission_tc: types.TheanoMatrix,
                                    resolve_nans: bool):
        """Generates symbolic tensors for log posterior and log data likelihood.
        
        Returns:
            log_posterior_probs, log_data_likelihood
        """

        def calculate_next_alpha(c_log_trans_mat: types.TheanoMatrix,
                                 c_log_emission_vec: types.TheanoVector,
                                 p_alpha_vec: types.TheanoVector):
            """Calculates the next entry on the forward table, alpha_{t}, from alpha_{t-1}.

            Args:
                c_log_trans_mat: a 2d tensor with rows and columns corresponding to log transition probability
                    from the previous state at position t-1 and to the next state at position t, respectively
                c_log_emission_vec: a 1d tensor representing the emission probability to each state at position t
                p_alpha_vec: a 1d tensor representing alpha_{t-1}

            Returns:
                symbolic 1d tensor of alpha_{t}
            """
            mu = tt.tile(p_alpha_vec, (num_states, 1)) + c_log_trans_mat.T
            n_alpha_vec = c_log_emission_vec + pm.math.logsumexp(mu, axis=1).dimshuffle(0)
            if resolve_nans:
                return tt.switch(tt.isnan(n_alpha_vec), -np.inf, n_alpha_vec)
            else:
                return n_alpha_vec

        def calculate_prev_beta(n_log_trans_mat: types.TheanoMatrix,
                                n_log_emission_vec: types.TheanoVector,
                                n_beta_vec: types.TheanoVector):
            """Calculates the previous entry on the backward table, beta_{t-1}, from beta_{t}.

            Args:
                n_log_trans_mat: a 2d tensor with rows and columns corresponding to log transition probability
                    from the previous state at position t-1 and to the next state at position t, respectively
                n_log_emission_vec: a 1d tensor representing the emission probability to each state at position t
                n_beta_vec: a 1d tensor representing beta_{t}

            Returns:
                symbolic 1d tensor of beta_{t-1}
            """
            nu = tt.tile(n_beta_vec + n_log_emission_vec, (num_states, 1)) + n_log_trans_mat
            p_beta_vec = pm.math.logsumexp(nu, axis=1).dimshuffle(0)
            if resolve_nans:
                return tt.switch(tt.isnan(p_beta_vec), -np.inf, p_beta_vec)
            else:
                return p_beta_vec

        # calculate thermal equivalent of various quantities
        inv_temperature = tt.inv(temperature)
        thermal_log_prior_c = inv_temperature * log_prior_c
        thermal_log_prior_c -= pm.math.logsumexp(thermal_log_prior_c)
        thermal_log_trans_tcc = inv_temperature * log_trans_tcc
        thermal_log_trans_tcc -= pm.math.logsumexp(thermal_log_trans_tcc, axis=-1)
        thermal_log_emission_tc = inv_temperature * log_emission_tc

        # first entry of the forward table
        alpha_first = thermal_log_prior_c + thermal_log_emission_tc[0, :]

        # the rest of the forward table
        alpha_seq, alpha_updates = th.scan(
            fn=calculate_next_alpha,
            sequences=[thermal_log_trans_tcc, thermal_log_emission_tc[1:, :]],
            outputs_info=[alpha_first])

        # concatenate with the first alpha
        alpha_t = tt.concatenate((alpha_first.dimshuffle('x', 0), alpha_seq))

        # last entry of the backward table (zero for all states)
        beta_last = tt.zeros_like(log_prior_c)

        # the rest of the backward table
        beta_seq, beta_updates = th.scan(
            fn=calculate_prev_beta,
            sequences=[thermal_log_trans_tcc, thermal_log_emission_tc[1:, :]],
            go_backwards=True,
            outputs_info=[beta_last])

        # concatenate with the last beta and reverse
        beta_t = tt.concatenate((beta_last.dimshuffle('x', 0), beta_seq))[::-1, :]

        # calculate normalized log posterior
        log_unnormalized_posterior_t = alpha_t + beta_t
        log_data_likelihood_t = pm.math.logsumexp(log_unnormalized_posterior_t, axis=1)
        log_posterior_t = log_unnormalized_posterior_t - log_data_likelihood_t

        return log_posterior_t, log_data_likelihood_t.dimshuffle(0), alpha_t, beta_t


class TheanoViterbi:
    """Implementation of Viterbi algorithm using `theano.scan`."""
    def __init__(self):
        self._viterbi_theano_func = self._get_compiled_viterbi_theano_func()

    def get_viterbi_path(self,
                         temperature: float,
                         log_prior_c: np.ndarray,
                         log_trans_tcc: np.ndarray,
                         log_emission_tc: np.ndarray):
        return self._viterbi_theano_func(temperature, log_prior_c, log_trans_tcc, log_emission_tc)

    @th.configparser.change_flags(compute_test_value="ignore")
    def _get_compiled_viterbi_theano_func(self):
        """todo.

        Args:

        Returns:
        """
        temperature = tt.scalar('temperature')
        log_prior_c = tt.vector('log_prior_c')
        log_trans_tcc = tt.tensor3('log_trans_tcc')
        log_emission_tc = tt.matrix('log_emission_tc')

        inputs = [temperature, log_prior_c, log_trans_tcc, log_emission_tc]
        viterbi_constrained_path_full_t = self._get_symbolic_log_viterbi_chain(
            temperature, log_prior_c, log_trans_tcc, log_emission_tc)

        return th.function(inputs=inputs, outputs=viterbi_constrained_path_full_t)

    @staticmethod
    def _get_symbolic_log_viterbi_chain(temperature: tt.scalar,
                                        log_prior_c: types.TheanoVector,
                                        log_trans_tcc: types.TheanoTensor3,
                                        log_emission_tc: types.TheanoMatrix):
        """Generates a symbolic 1d integer tensor representing the most likely chain of hidden states
        (Viterbi algorithm).

        Returns:
            symbolic 1d integer tensor representing the most likely chain of hidden states
        """

        def calculate_next_omega_psi(p_log_trans_ab, c_log_emission_b, p_omega_a):
            """Calculates the updated data log likelihood and the next best state by taking into the
            next link in the Markov chain.

            Args:
                p_log_trans_ab: log transition matrix from `a` to `b`
                c_log_emission_b: log emission probabilities at the current position
                p_omega_a: previous data log likelihood assuming the last state is `a`

            Returns:
                updated data log likelihood for all possible terminal states,
                previous best state conditioned on each of the states for the current position
            """
            tau_ab = p_log_trans_ab + p_omega_a.dimshuffle(0, 'x')
            max_tau_b, best_state_b = tt.max_and_argmax(tau_ab, axis=0)
            n_omega_b = c_log_emission_b + max_tau_b
            return n_omega_b, best_state_b

        def calculate_previous_best_state(c_psi_c, c_best_state):
            """Calculates the previous best state (trivially) from the current best state and
            the Viterbi table (psi).

            Args:
                c_psi_c: Viterbi table (previous best state conditioned on the current best state)
                c_best_state: best state in the current position

            Returns:
                previous best state
            """
            return c_psi_c[c_best_state]

        # calculate thermal equivalent of various quantities
        inv_temperature = tt.inv(temperature)
        thermal_log_prior_c = inv_temperature * log_prior_c
        thermal_log_prior_c -= pm.math.logsumexp(thermal_log_prior_c)
        thermal_log_trans_tcc = inv_temperature * log_trans_tcc
        thermal_log_trans_tcc -= pm.math.logsumexp(thermal_log_trans_tcc, axis=-1)
        thermal_log_emission_tc = inv_temperature * log_emission_tc

        # state log likelihood for the first position
        omega_first_a = thermal_log_emission_tc[0, :] + thermal_log_prior_c

        # calculate the the log likelihood of partial Viterbi paths (omega_tc) and the Viterbi table (psi_tc)
        omega_psi_list, _ = th.scan(
            fn=calculate_next_omega_psi,
            sequences=[thermal_log_trans_tcc, thermal_log_emission_tc[1:, :]],
            outputs_info=[omega_first_a, None])
        omega_tc = omega_psi_list[0]
        psi_tc = omega_psi_list[1]

        # obtain the Viterbi chain from omega_tc and psi_tc
        last_best_state = tt.argmax(omega_tc[-1, :])
        viterbi_constrained_path_except_for_last_t, _ = th.scan(
            fn=calculate_previous_best_state,
            sequences=[psi_tc],
            outputs_info=[last_best_state],
            go_backwards=True)

        viterbi_constrained_path_full_t = tt.concatenate(
            [tt.stack(last_best_state),
             viterbi_constrained_path_except_for_last_t])[::-1]

        return viterbi_constrained_path_full_t


class HMMSegmentationQualityCalculator:
    """Calculates quality metrics for hidden state segments for a given HMM.

    Note:
        The initializer requires the emission and transition probabilities, as well as the forward
        and backward tables and the log posterior probability.
    """

    # 10 / ln(10)
    INV_LN_10_TIMES_10 = 4.342944819032518

    # ln(1/2)
    LN_HALF = -0.6931471805599453

    def __init__(self,
                 log_emission_tc: np.ndarray,
                 log_trans_tcc: np.ndarray,
                 alpha_tc: np.ndarray,
                 beta_tc: np.ndarray,
                 log_posterior_prob_tc: np.ndarray,
                 log_data_likelihood: float):
        assert isinstance(log_emission_tc, np.ndarray)
        assert log_emission_tc.ndim == 2
        self.num_sites, self.num_states = log_emission_tc.shape

        assert isinstance(log_trans_tcc, np.ndarray)
        assert log_trans_tcc.shape == (self.num_sites - 1, self.num_states, self.num_states)

        assert isinstance(alpha_tc, np.ndarray)
        assert alpha_tc.shape == (self.num_sites, self.num_states)

        assert isinstance(beta_tc, np.ndarray)
        assert beta_tc.shape == (self.num_sites, self.num_states)

        assert isinstance(log_posterior_prob_tc, np.ndarray)
        assert log_posterior_prob_tc.shape == (self.num_sites, self.num_states)

        self.log_emission_tc = log_emission_tc
        self.log_trans_tcc = log_trans_tcc
        self.alpha_tc = alpha_tc
        self.beta_tc = beta_tc
        self.log_posterior_prob_tc = log_posterior_prob_tc
        self.log_data_likelihood = log_data_likelihood

        self._constrained_path_log_likelihood_theano_func =\
            self._get_compiled_constrained_path_log_likelihood_theano_func()

        self.all_states_set = set(range(self.num_states))

    @staticmethod
    def _get_symbolic_constrained_path_log_likelihood(alpha_first_c: tt.dvector,
                                                      beta_last_c: tt.dvector,
                                                      log_emission_tc: tt.dmatrix,
                                                      log_trans_tcc: tt.dtensor3) -> tt.dscalar:

        def update_alpha(c_log_emission_c: tt.dvector,
                         c_log_trans_cc: tt.dmatrix,
                         p_alpha_c: tt.dvector):
            return c_log_emission_c + pm.math.logsumexp(p_alpha_c.dimshuffle(0, 'x') +
                                                        c_log_trans_cc, axis=0).dimshuffle(1)

        alpha_seg_iters, _ = th.scan(
            fn=update_alpha,
            sequences=[log_emission_tc, log_trans_tcc],
            outputs_info=[alpha_first_c])
        alpha_seg_end_c = alpha_seg_iters[-1, :]

        return pm.math.logsumexp(alpha_seg_end_c + beta_last_c)

    @staticmethod
    @th.configparser.change_flags(compute_test_value="ignore")
    def _get_compiled_constrained_path_log_likelihood_theano_func():
        alpha_first_c = tt.dvector('alpha_first_c')
        beta_last_c = tt.dvector('beta_last_c')
        log_emission_tc = tt.dmatrix('log_emission_tc')
        log_trans_tcc = tt.dtensor3('log_trans_tcc')

        inputs = [alpha_first_c, beta_last_c, log_emission_tc, log_trans_tcc]
        output = HMMSegmentationQualityCalculator._get_symbolic_constrained_path_log_likelihood(
            alpha_first_c, beta_last_c, log_emission_tc, log_trans_tcc)

        return th.function(inputs=inputs, outputs=output)

    def get_log_constrained_posterior_prob(self,
                                           start_index: int, end_index: int,
                                           allowed_states: Set[int]) -> float:
        """Calculates the constrained log posterior probability for contiguous set of sites in a Markov chain.
        At each site, only a subset of all states (as set by `allowed_states`) are allowed and the other states
        are strictly avoided.

        Args:
            start_index: first site index (inclusive)
            end_index: last site index (inclusive)
            allowed_states: the list of allowed states in the segment

        Returns:
            log constrained posterior probability (float)
        """
        assert start_index >= 0
        assert end_index < self.num_sites
        assert end_index >= start_index
        assert all(isinstance(item, int) and 0 <= item < self.num_states for item in allowed_states),\
            "The set of allowed states must be integers and in range [0, {0}]".format(self.num_states - 1)
        allowed_states_list = sorted(allowed_states)

        constrained_alpha_first_c = self.alpha_tc[start_index, allowed_states_list]
        constrained_beta_last_c = self.beta_tc[end_index, allowed_states_list]
        if end_index == start_index:  # single-site segment
            log_constrained_data_likelihood: float = logsumexp(constrained_alpha_first_c + constrained_beta_last_c)
            return log_constrained_data_likelihood - self.log_data_likelihood
        else:
            # calculate the required slices of the log emission and log transition representing
            # paths that only contain the allowed states
            constrained_log_emission_tc =\
                self.log_emission_tc[(start_index + 1):(end_index + 1), allowed_states_list]
            constrained_log_trans_tcc =\
                self.log_trans_tcc[start_index:end_index, allowed_states_list, :][:, :, allowed_states_list]
            return np.asscalar(self._constrained_path_log_likelihood_theano_func(
                constrained_alpha_first_c, constrained_beta_last_c,
                constrained_log_emission_tc, constrained_log_trans_tcc) - self.log_data_likelihood)

    def get_segment_some_quality(self, start_index: int, end_index: int, call_state: int) -> float:
        """Calculates the phred-scaled posterior probability that one or more ("some") sites in a segment have
        the same hidden state ("call").

        Args:
            start_index: first site index (inclusive)
            end_index: last site index (inclusive)
            call_state: segment call state index

        Returns:
            a phred-scaled probability
        """
        assert call_state in self.all_states_set
        other_states = self.all_states_set.copy()
        other_states.remove(call_state)
        all_other_states_log_prob = self.get_log_constrained_posterior_prob(start_index, end_index, other_states)
        return self.log_prob_to_phred(all_other_states_log_prob, complement=False)

    def get_segment_exact_quality(self, start_index: int, end_index: int, call_state: int) -> float:
        """Calculates the phred-scaled posterior probability that "all" sites in a segment have the same
        hidden state ("call").

        Args:
            start_index: first site index (inclusive)
            end_index: last site index (inclusive)
            call_state: segment call state index

        Returns:
            a phred-scaled probability
        """
        assert call_state in self.all_states_set
        all_called_state_log_prob = self.get_log_constrained_posterior_prob(start_index, end_index, {call_state})
        return self.log_prob_to_phred(all_called_state_log_prob, complement=True)

    def get_segment_start_quality(self, start_index: int, call_state: int) -> float:
        """Calculates the phred-scaled posterior probability that a site marks the start of a segment.

        Args:
            start_index: left breakpoint index of a segment
            call_state: segment call state index

        Returns:
            a phred-scaled probability
        """
        assert 0 <= start_index < self.num_sites
        if start_index == 0:
            log_prob = self.log_posterior_prob_tc[0, call_state]
        else:
            # calculate the probability of all paths that start from other states and end up with the called state
            all_other_states = self.all_states_set.copy()
            all_other_states.remove(call_state)
            all_other_states_list = list(all_other_states)
            prev_alpha_c = self.alpha_tc[start_index - 1, all_other_states_list]
            current_beta = self.beta_tc[start_index, call_state]
            current_log_emission = self.log_emission_tc[start_index, call_state]
            log_trans_c = self.log_trans_tcc[start_index - 1, all_other_states_list, call_state]
            log_breakpoint_likelihood = logsumexp(prev_alpha_c + log_trans_c + current_log_emission) + current_beta
            log_prob = log_breakpoint_likelihood - self.log_data_likelihood

        return self.log_prob_to_phred(log_prob, complement=True)

    def get_segment_end_quality(self, end_index: int, call_state: int) -> float:
        """Calculates the phred-scaled posterior probability that a site marks the end of a segment.

        Args:
            end_index: right breakpoint index of a segment
            call_state: segment call state index

        Returns:

        """
        assert 0 < end_index < self.num_sites
        if end_index == self.num_sites - 1:
            log_prob = self.log_posterior_prob_tc[self.num_sites - 1, call_state]
        else:
            # calculate the probability of all paths that start from call state and end up with other states
            all_other_states = self.all_states_set.copy()
            all_other_states.remove(call_state)
            all_other_states_list = list(all_other_states)
            current_alpha = self.alpha_tc[end_index, call_state]
            next_beta_c = self.beta_tc[end_index + 1, all_other_states_list]
            next_log_emission_c = self.log_emission_tc[end_index + 1, all_other_states_list]
            log_trans_c = self.log_trans_tcc[end_index, call_state, all_other_states_list]
            log_breakpoint_likelihood = logsumexp(current_alpha + log_trans_c + next_log_emission_c + next_beta_c)
            log_prob = log_breakpoint_likelihood - self.log_data_likelihood

        return self.log_prob_to_phred(log_prob, complement=True)

    @staticmethod
    def log_prob_to_phred(log_prob: float, complement: bool = False) -> float:
        """Converts probabilities in natural log scale to phred scale.

        Args:
            log_prob: a probability in the natural log scale
            complement: invert the probability

        Returns:
            phred-scaled probability
        """
        final_log_prob = log_prob if not complement else HMMSegmentationQualityCalculator.log_prob_complement(log_prob)
        return -final_log_prob * HMMSegmentationQualityCalculator.INV_LN_10_TIMES_10

    @staticmethod
    def log_prob_complement(log_prob: float) -> float:
        """Calculates the complement of a probability in the natural log scale.

        Args:
            log_prob: a probability in the natural log scale

        Returns:
            complement of the the probability in the natural log scale
        """
        log_prob_zero_capped = min(0., log_prob)
        if log_prob_zero_capped >= HMMSegmentationQualityCalculator.LN_HALF:
            return np.log(-np.expm1(log_prob_zero_capped))
        else:
            return np.log1p(-np.exp(log_prob_zero_capped))
