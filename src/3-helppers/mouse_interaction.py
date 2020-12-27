# force based

    # constraint mouse selection
    if mouse["selected"]:
        f_ex[mouse["selected"]] = 750 * (mouse["x"] - x[mouse["selected"]])

    # free mouse mouse selection
    if mouse["selected"]:
        f_ex[mouse["selected"]] = [0, 0]

# position based

    # constraint mouse selection
    if mouse["selected"]:
        free = np.setdiff1d(free, mouse["selected"])
        x[mouse["selected"]] = mouse["x"]

    # free mouse mouse selection
    if mouse["selected"]:
        free = np.setdiff1d(np.arange(mesh.ndof), cons)
