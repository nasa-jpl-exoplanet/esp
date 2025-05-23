'''cerberus bot ds'''

# pylint: disable=abstract-method,arguments-differ

# -- IMPORTS -- ------------------------------------------------------
import numpy as np

import pytensor.graph as tnsrgraph
import pytensor.tensor as tnsr


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
    ForwardModel = orbital(*newnodes)
    # Norm = np.log(np.sqrt(2e0 * np.pi)) - np.log(ctxt.mcmcsig)
    Norm = np.log(2e0 * np.pi * ctxt.mcmcsig)
    out = -(((ctxt.mcmcdat - ForwardModel) / ctxt.mcmcsig) ** 2) / 2e0 - Norm

    return out


class TensorShell(LogLikelihood, tnsrgraph.Op):
    '''
    GMR: Tensor Shell for custom models
    Do not touch the name of the methods
    '''

    def make_node(self, nodes) -> tnsrgraph.Apply:
        inputs = [tnsr.as_tensor(n) for n in nodes]
        outputs = [tnsr.vector()]
        return tnsrgraph.Apply(self, inputs, outputs)

    def perform(
        self,
        node: tnsrgraph.Apply,
        inputs: list[np.ndarray],
        output_storage: list[list[None]],
    ) -> None:
        output_storage[0][0] = np.asarray(LogLikelihood(inputs))
        return

    pass
