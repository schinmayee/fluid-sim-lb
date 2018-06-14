class Configuration(object):
    # Instance group
    group = None
    # Zone
    zone = None
    # Worker threads
    worker_threads = 8
    # Initial placement, before the simulation even starts
    placement = "application"
    # Affinity parameters
    ax = 4
    ay = 2
    az = 2
    # Simulation scale and number of partitions
    ngrid = 1024
    x = 8
    y = 8
    z = 8
    horizon = 30  # controls frequency of load balancing
    current_load_factor = 0.0
    use_geometric = False
    # Simulation parameters
    init = 0  # initial configuration
    cfl = 8  # cfl number
    frames = 200  # number of frames
    p = 1 # number of processors per worker, igngore this for now
    flip = 0.4  # flip factor
    frame_rate = 30  # frame rate
    solver_iterations = 2000  # maximum number of solver iterations
    # Coarse simulation scaling factor and partitions if applicable
    coarse_scale = 8
    cx = 2
    cy = 2
    cz = 2
    window = 1  # controls how far ahead coarse simulation runs
    init_file = ""

"""
Parse and return experiment configuration.
"""
def ParseConfiguration(config_file):
    config = Configuration()
    with open(config_file) as data:
        for line in data:
            if len(line) == 0:
                continue
            if line[0] == "#":
                # Comments
                continue
            print(line)
            tokens = line.split(':')
            key = tokens[0].strip()
            val = tokens[1].strip()
            if key == "group":
                config.group = val
            elif key == "zone":
                config.zone = val
            elif key == "worker_threads":
                config.worker_threads = int(val)
            elif key == "placement":
                config.placement = val
            elif key == "ax":
                config.ax = int(val)
            elif key == "ay":
                config.ay = int(val)
            elif key == "az":
                config.az = int(val)
            elif key == "ngrid":
                config.ngrid = int(val)
            elif key == "x":
                config.x = int(val)
            elif key == "y":
                config.y = int(val)
            elif key == "z":
                config.z = int(val)
            elif key == "horizon":
                config.horizon = int(val)
            elif key == "current_load_factor":
                config.current_load_factor = float(val)
            elif key == "coarse_scale":
                config.coarse_scale = int(val)
            elif key == "cx":
                config.cx = int(val)
            elif key == "cy":
                config.cy = int(val)
            elif key == "cz":
                config.cz = int(val)
            elif key == "window":
                config.window = int(val)
            elif key == "init":
                config.init = int(val)
            elif key == "cfl":
                config.cfl = int(val)
            elif key == "frames":
                config.frames = int(val)
            elif key == "p":
                config.p = int(val)
            elif key == "flip":
                config.flip = float(val)
            elif key == "frame_rate":
                config.frame_rate = int(val)
            elif key == "solver_iterations":
                config.solver_iterations = int(val)
            elif key == "use_geometric":
                config.use_geometric = bool(val)
            elif key == "init_file":
                config.init_file = val

    assert(config.group is not None)
    assert(config.zone is not None)
    attrs = vars(config)
    print '\n'.join("%s: %s" % item for item in attrs.items())
    return config
