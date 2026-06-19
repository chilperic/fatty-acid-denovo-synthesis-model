import numpy as np
from fadns_model.simulation import simulate_time_course


def toy_rhs(t, y, params):
    k = params.get("k", 0.2)
    return np.array([-k * y[0], k * y[0] - 0.1 * y[1]])

result = simulate_time_course(toy_rhs, y0=[1.0, 0.0], params={"k": 0.2}, t_span=(0, 10), n_points=50)
print("success:", result.success)
print("final_state:", result.y[:, -1])
