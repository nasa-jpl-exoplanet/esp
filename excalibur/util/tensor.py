'''util tensor ds'''

# Heritage code shame:
# pylint: disable=invalid-name

import numpy as np

import pytensor.graph as tnsrgraph
import pytensor.tensor as tnsr

from excalibur.cerberus.fmcontext import ctxtinit


# this doesn't change results at all; just needed to avoid undefined-variable pylint
ctxt = ctxtinit()


class TensorShell(tnsrgraph.Op):
    '''
    GMR: Tensor Shell for custom models
    Do not touch the name of the methods
    GB: R_op and grad definitions added to avoid abstract-method pylint error
    '''

    def make_node(self, *nodes) -> tnsrgraph.Apply:
        inputs = [tnsr.as_tensor(n) for n in nodes[0]]
        outputs = [tnsr.vector()]
        return tnsrgraph.Apply(self, inputs, outputs)

    def R_op(self, *_args, **_keywords):
        raise NotImplementedError('not expecting this method to be used')

    def grad(self, *_args, **_keywords):
        raise NotImplementedError('not expecting this method to be used')

    def perform(
        self,
        node: tnsrgraph.Apply,
        inputs: list[np.ndarray],
        output_storage: list[list[None]],
    ) -> None:
        output_storage[0][0] = np.asarray(LogLikelihood(inputs))
        return

    pass


# GMR: Gregoire s legacy
def LogLikelihood(inputs):
    '''
    GMR: User defined loglikelihood
    We stick to the proper definition of it
    '''
    newnodes = []
    newindex = 0
    for ns in ctxt.nodeshape:
        if ns > 1:
            newnodes.append(inputs[newindex : newindex + ns])
            pass
        else:
            newnodes.append(inputs[newindex])
            pass
        newindex += ns
        pass
    # ForwardModel = orbital(*newnodes)
    # ForwardModel = clearfmcerberus(*newnodes)
    # ForwardModel = cloudyfmcerberus(*newnodes)
    ForwardModel = ctxt.forwardmodel(*newnodes)

    out = -(((ctxt.mcmcdat - ForwardModel) / ctxt.mcmcsig) ** 2) / 2e0

    # this is a useful check; chi2_red should decrease toward ~1 (for simulated data)
    # print('  chi2_reduced for this model:', -2 * np.sum(out) / len(out))

    # normalize the log(Likelihood); as a constant, it shouldn't have any effect
    Norm = np.log(2e0 * np.pi * ctxt.mcmcsig)
    out -= Norm

    return out
