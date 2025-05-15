'''cerberus bot ds'''

# -- IMPORTS -- ------------------------------------------------------
import pytensor.graph as tnsrgraph
import pytensor.tensor as tnsr

class TensorShell(tnsrgraph.Op):
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
