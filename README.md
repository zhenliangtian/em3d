This program simulates the evolution of the Earth-Moon-Sun system. It integrates the gravitational forces between point masses, as well as the free rigid body rotations of Earth and Moon, the spin-orbit interactions between the Earth's figure and the other two bodies, the interaction between the Moon's figure with the Earth, and the tidal interactions between the Earth and Moon. Earth's figure varies depending on the rotation rate of the Earth. The constant-Q Darwin-Kaula model is used for tides. The symplectic mapping method is used for the N-body integration.


The procedures to run the program:

1. Download the program files (main.c, *.h) and save them in a directory: ./;

2. Compile in shell: gcc main.c -lm;

3. Create a data directory: mkdir ./data/;

4. Run the executable, with a particular run id:
    ./a.out ANYID;

5. In the console, the program will ask the used to set system parameters (k2 and Q for the Earth and Moon), the initial condition (Earth-Moon distance, lunar eccentricity, etc.), and computation control parameters (the range of lunar semi-major axis, outside which the program will self-terminate; and parameters that control how frequent the program should write the system state into the data file, and the maximum number of steps of integration).

6. Once the initial condition and parameter values are set, the program starts the integration, and will continue to write system states into the data file until either the lunar orbit's semi-major axis exceeds the given bounds, or the maximum number of steps of integration is reached.


After the integration terminates, there will be four files in ./data/:
    meta_ANYID.txt     (metadata recording the time of beginning and termination, and other details of the integration),
    ic_ANYID.txt       (record of the initial condition and the tidal parameters),
    evol_ANYID.out     (binary data containing system states with the specialized interval),
    resuinfo_ANYID.out (binary data that is helpful in case the user would like to extend or resume the run later).


This program is made public on Apr 29, 2020.
