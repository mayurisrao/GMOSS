from pygmoss import utils
import warnings

warnings.filterwarnings("ignore")


def test_concave_func():
    assert (
        utils.concave_func(
            1.7,
            -1.7673603339907888,
            6.209947892799811,
            2.8215667328790097,
            2.6599927540080444,
            0.16663361137887173,
            9999.790968653713,
            0.00020113422862317362,
        )
        == 1.0218574884906784
    )


def test_convex_func():
    assert (
        utils.convex_func(
            1.7,
            8.317511113295733e-07,
            2.51968636124662,
            2.707635701950907,
            247386337.5370596,
            7338.466282609452,
            8000.0,
            0.001,
        )
        == 2.7367684788205375
    )
