# Forward Euler x
    y = x + v * dt + M_inv.dot(f_ex) * dt**2
    x[free] = (y + M_inv.dot(compute_f(x)) * dt**2)[free]

    # Forward Euler as dx
    y = x + v * dt + f_ex * dt**2
    dx = (y - x + M_inv.dot(compute_f(x)) * dt**2)[free]
    x[free] += dx

    # Backward Euler as x
    y = x + v * dt + M_inv.dot(f_ex) * dt**2
    a = M - compute_k(x) * dt**2
    b = M[free].T[free].T.dot(y[free]) + compute_f(x)[free] * dt**2 - compute_k(x)[free].T[free].T.dot(x[free]) * dt**2
    x[free] = np.linalg.solve(a[free].T[free].T, b)

    x[free] = np.linalg.solve(a, b)[free]  # is not accurate

    # Backward Euler as dx
    y = x + v * dt + M_inv.dot(f_ex) * dt**2
    a = M - compute_k(x) * dt**2
    b = M.dot(y - x) + compute_f(x) * dt**2
    dx = np.linalg.solve(a[free].T[free].T, b[free])
    x[free] += dx

    dx = np.linalg.solve(a, b)  # is not accurate
    x[free] += dx[free]  # is not accurate

    # Static as dx
    a = -compute_k(x)
    b = compute_f(x) + f_ex
    dx = np.linalg.solve(a[free].T[free].T, b[free])
    x[free] += dx

    # Newmark as x
    a = np.identity(mesh.ndof) - compute_k(x) * dt**2 / 4
    y = x + v * dt + compute_f(x) * dt**2 / 4 + f_ex * dt**2 / 2
    b = y[free] + compute_f(x)[free] * dt**2 / 4 - compute_k(x)[free].T[free].T.dot(x[free]) * dt**2 / 4
    x[free] = np.linalg.solve(a[free].T[free].T, b)
    v = (x - xn - v * dt / 2) * (2 / dt)
