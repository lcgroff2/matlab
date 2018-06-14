%LatGasDiff.m: Nanoparticle swelling/energy transfer/exciton diffusion
%simulation based on lattice gas model -LCG 02/03/2014

%Summary: Random walk model simulating exciton motion in a cubic
%nanoparticle. Accounts for changes in chromophore density as
%solvent quality shifts, where excitons hop from one
%available chromophore to the nearest neighbor with each timestep.
%Excitons may either hop between neighbor chromophores, decay
%radiatively, or undergo energy transfer to a quencher. Here, a
%quencher is defined as either a redshifted aggregate chromophore
%or a nonfluorescent chromophore.

%Adjustable Parameters:

%Nanoparticle Radius (nm):
param.nprad = 5;
%Exciton diffusion length L_D (nm):
param.ld = 3;
%Polymer-polymer Förster Radius (nm)
param.r0 = 4; %Assume perylene red-like quenchers for now.
%Cubic lattice spacing (nm):
param.latspc = 0.1;
%Number of chromophores:
param.chrom = round(2*param.nprad./param.latspc)^3;
%Number of quenchers:
param.quench = 0;
%polymer lifetime in good solvent (ps):
param.tau0 = 3000;
%Number of excitons to populate in nanoparticle:
param.nex = 1000;
%Number of averages (individual nanoparticles) to compute:
param.navg = 1;
%Time step (ps)
param.dt = 0.1;

%Diffusion parameters:
D = param.ld^2/(6*param.tau0); %Diffusion Constant (3D)
sig = sqrt(2*D*param.dt); %Theoretical RMS Displacement
steps = round(6*param.tau0/param.dt); %Number of steps in each individual walk
prob = sig^2/param.latspc^2; %Exciton hopping probability
disp(param)
fprintf(1, 'Exciton Hopping Probability = %1.3f \n',prob)

%A Few Checks:
if param.chrom > (2*param.nprad/param.latspc)^3
    error('Too many chromophores!')
end

if (param.chrom+param.quench) > (2*param.nprad/param.latspc)^3
    error('Too many chromophores + Quenchers!')
end

results.sumpop = zeros([1 steps]);
for navg = 1:param.navg,
    % Build Lattice as 3D empty matrix, lattice points are denoted by indeces
    % on vectors:
    np.latpts = zeros(round(2*param.nprad/param.latspc+2), round(2*param.nprad/param.latspc+2), round(2*param.nprad/param.latspc+2));
    %for 10 nm particle diameter, 0.1 nm spacing =  100x100x100 matrix of 10^6
    %lattice points. Add two to each dimension to make it 102x102x102 matrix.
    %Randomly generate chromophores/quenchers within indeces 2:101 to confine
    %the 100^3 matrix within a wall of zeros to provide a boundary for the
    %excitons in case they propagate onto the faces of the nanoparticle.

    %Pre-allocate empty matrices for keeping track of chromophore/quencher
    %locations. Will need to store a copy of each chromophore position, and(
    %each quencher position separately, then concatenate the vectors together
    %to make the combined list of possible positions to initially place
    %excitons.
    np.cposx = zeros(1, length(param.chrom));
    np.cposy = np.cposx;
    np.cposz = np.cposx;
    np.qposx = zeros(1, length(param.quench));
    np.qposy = np.qposx;
    np.qposz = np.qposx;



    %Insert Chromophores:
    %Lattice points are integers in a 1x100 vector for each direction. Use
    %randint to generate the vector index positions where a 1 is placed to
    %indicate a chromophore.

    rand('state', sum(100*clock)) %Randomize the RNG
    cnt = 1; %dummy variable
    while cnt <= param.chrom,
        n = round(2*param.nprad/param.latspc);
        np.xtry = ceil(n.*rand)+1; %Generate random integers from 1:n, +1 so they stay in the interval 2:n+1;
        np.ytry = ceil(n.*rand)+1;
        np.ztry = ceil(n.*rand)+1;
        if np.latpts(np.xtry, np.ytry, np.ztry) == 0
            np.latpts(np.xtry, np.ytry, np.ztry) = 1;  %Place chromophore in NP
            %Save a copy of chromophore position
            np.cposx(cnt) = np.xtry;
            np.cposy(cnt) = np.ytry;
            np.cposz(cnt) = np.ztry;
            cnt = cnt+1;
        end
    end

    %Sanity check to make sure chromophores were placed properly
    chromsum = sum(np.latpts); %3d matrix to 2d matrix
    chromsum = sum(chromsum);  %2d matrix to 1d vector
    chromsum = sum(chromsum);  %1d vector to scalar
    if chromsum == param.chrom
        disp('Chromophores placed successfully.')
    else
        disp('Something went wrong.')
    end

    %Insert Quenchers:
    cnt2 = 1; %dummy variable
    if param.quench >= 1
        while cnt2 <= param.quench
            n = round(2*param.nprad/param.latspc);
            q.xtry = ceil(n.*rand)+1; %Generate random integers from 1:n, +1 so they stay in the interval 2:n+1;
            q.ytry = ceil(n.*rand)+1;
            q.ztry = ceil(n.*rand)+1;
            %If no chromophore or quencher exists at this point put quencher there.
            if np.latpts(q.xtry, q.ytry, q.ztry) == 0
                np.latpts(q.xtry, q.ytry, q.ztry) = -1;
                %Store a copy of each quencher position
                np.qposx(cnt2) = q.xtry;
                np.qposy(cnt2) = q.ytry;
                np.qposz(cnt2) = q.ztry;
                cnt2 = cnt2+1;
            end
        end
    end
    qsum = sum(np.latpts);
    qsum = sum(qsum);
    qsum = sum(qsum);
    if qsum == param.chrom-param.quench
        disp('quenchers placed successfully.')
    else
        disp('something went wrong.')
    end
    %Concatenate position vectors together for chromophores and quenchers:
    if param.quench >= 1
        np.cqposx = horzcat(np.cposx, np.qposx);
        np.cqposy = horzcat(np.cposy, np.qposy);
        np.cqposz = horzcat(np.cposz, np.qposz);
    else
        np.cqposx = np.cposx;
        np.cqposy = np.cposy;
        np.cqposz = np.cposz;
    end
    %Place excitons on chromophores/quenchers:


    decayed = zeros([1 param.nex]);
    alive = decayed+1;
    quenched = decayed;
    decaytime = decayed;

    cnt3 = 1;
    while cnt3 <= param.nex

        %Physical picture:

        %Treat excitons as independent, non-interacting units, so that we can have
        %more excitons than chromophores+quenchers for a given nanoparticle. The
        %idea is that chromophores = 1's in np matrix, add 1 to each chromophore
        %where we place an exciton, regardless of whether or not an exciton is
        %already there. As they move/decay, we stop the downward count on a lattice
        %point once it returns to 1. If an exciton is placed on a quencher,
        %subtract by 1 for each exciton placed on that quencher, and the
        %upward count stops at -1 as they decay. At time zero, remove these
        %quenched excitons from the population, place in a "quenched" list
        %and proceed.

        %Pick a random position index
        n2 = length(np.cqposx);
        np.extry = ceil(n2*rand);

        if np.latpts(np.cqposx(np.extry), np.cqposy(np.extry), np.cqposz(np.extry)) >= 1
            np.latpts(np.cqposx(np.extry), np.cqposy(np.extry), np.cqposz(np.extry)) = np.latpts(np.cqposx(np.extry), np.cqposy(np.extry), np.cqposz(np.extry))+1;
            np.exposx(cnt3) = np.cqposx(np.extry);
            np.exposy(cnt3) = np.cqposy(np.extry);
            np.exposz(cnt3) = np.cqposz(np.extry);
            cnt3 = cnt3+1;
        elseif np.latpts(np.cqposx(np.extry), np.cqposy(np.extry), np.cqposz(np.extry)) <= -1
            np.latpts(np.cqposx(np.extry), np.cqposy(np.extry), np.cqposz(np.extry)) = np.latpts(np.cqposx(np.extry), np.cqposy(np.extry), np.cqposz(np.extry))-1;
            np.exposx(cnt3) = np.cqposx(np.extry);
            np.exposy(cnt3) = np.cqposy(np.extry);
            np.exposz(cnt3) = np.cqposz(np.extry);
            cnt3 = cnt3+1;
        end
    end

    %Track population vs. time
    results.pop = zeros([1 steps]);
    results.time = results.pop;

    for curstep = 1:steps,
        t = param.dt*(curstep-1); %Time tracking

        results.pop(curstep) = sum(alive);
        results.time(curstep) = t;

        %Check for excitons placed on quenchers. Subtract them from the
        %population at time zero. (Static Quenching)
        if t == 0;
            cnt4 = 1:1:param.nex;
            for cnt4 = 1:param.nex,
                if np.latpts(np.exposx(cnt4), np.exposy(cnt4), np.exposz(cnt4)) < -1
                    np.latpts(np.exposx(cnt4), np.exposy(cnt4), np.exposz(cnt4)) = np.latpts(np.exposx(cnt4), np.exposy(cnt4), np.exposz(cnt4))+1;
                    quenched(cnt4) = 1;
                    alive(cnt4) = 0;
                    decaytime(cnt4) = t;
                end
            end
        end

        %See if all excitons have decayed. Exit loop if true:
        nlive = sum(alive);
        aliveidx = find(alive);
        if nlive < 1
            break
        end

        %Some excitons decay to ground state before hopping:
        justdecayed = alive .* (exp(-param.dt/param.tau0) < rand([1 param.nex]));
        alive = alive - justdecayed;
%         check = 1:1:param.nex;
%         for check = 1:param.nex,
%             if alive(check) < 0
%                 error('Line 229: alive vector has gone into negative numbers!')
%             end
%         end
        decayed = decayed + justdecayed;
        %Track decay times:
        decayidx = find(justdecayed==1);
        decaytime(decayidx) = t;

        %Some excitons undergo FRET before hopping:
        if param.quench >= 1
            aliveidx = find(alive);
            nlive = sum(alive);
            %Calculate quencher radial distances relative to exciton:
            rsq = zeros([1 nlive]);
            rates = zeros([1 nlive]);
            idx = 1:length(np.qposx);
            FRET = param.r0^6/param.tau0;
            for idx = 1:length(np.qposx),
                dx = (np.qposx(idx)-np.exposx(aliveidx)).*param.latspc;
                dy = (np.qposy(idx)-np.exposy(aliveidx)).*param.latspc;
                dz = (np.qposz(idx)-np.exposz(aliveidx)).*param.latspc;
                rsq = dx.*dx + dy.*dy + dz.*dz;
                rates = rates + FRET./(rsq.*rsq.*rsq);
            end
            %calculate FRET probability prior to hopping:
            justquenched = exp(-rates.*param.dt) < rand([1 nlive]);
            alive(aliveidx) = alive(aliveidx)-justquenched;
            quenched(aliveidx) = quenched(aliveidx)+justquenched;
            quenchidx = find(justquenched == 1);
            decaytime(quenchidx) = t;

        end
        if t == 0,
            pos0x = np.exposx;
            pos0y = np.exposy;
            pos0z = np.exposz;
        end
        %Assign a hop direction along each axis for all excitons, then move them.
        i = 1:1:3;
        for i = 1:3 %iterate 1x for each axis

            cnt5 = 1:1:param.nex;
            probhop = rand([1 param.nex]) < prob; %Generate random number, compare against hopping probability to see which excitons move.
            hopdir = floor(rand([1 param.nex])*2)*2-1; %Generate random direction for excitons to move. +1 or -1
            aliveidx = find(alive);
            exhop = probhop(aliveidx).*hopdir(aliveidx); %From the pool of alive excitons, determine whether they move or stay put.
            for cnt5 = 1:length(exhop),
                if exhop(cnt5) == 1
                    %Exciton hops forward, if possible.
                    switch i
                        case 1 %X axis
                            if np.latpts(np.exposx(aliveidx(cnt5))+1, np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5))) ~= 0
                                if np.latpts(np.exposx(aliveidx(cnt5))+1, np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5))) >= 1 %If a chromophore is present
                                    %Subtract exciton off of current position
                                    np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5))) = np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5)))-1;
                                    np.exposx(aliveidx(cnt5)) = np.exposx(aliveidx(cnt5))+1; %move exciton one step forward
                                    %Add exciton to forward position
                                    np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5))) = np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5)))+1;
                                elseif np.latpts(np.exposx(aliveidx(cnt5))+1, np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5))) == -1 %hopping to a quencher
                                    %Subtract exciton from current position
                                    np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5))) = np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5)))-1;
                                    quenched(aliveidx(cnt5)) = 1;
                                    alive(aliveidx(cnt5)) = 0;
                                    decaytime(aliveidx(cnt5)) = t;
                                end
                            end %If exciton neighbor = 0, cannot move forward, exciton stays put.
                        case 2 %Y axis
                            if np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5))+1, np.exposz(aliveidx(cnt5))) ~= 0
                                if np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5))+1, np.exposz(aliveidx(cnt5))) >= 1 %If a chromophore is present
                                    %Subtract exciton off of current position
                                    np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5))) = np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5)))-1;
                                    np.exposy(aliveidx(cnt5)) = np.exposy(aliveidx(cnt5))+1; %move exciton one step forward
                                    %Add exciton to forward position
                                    np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5))) = np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5)))+1;
                                elseif np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5))+1, np.exposz(aliveidx(cnt5))) == -1 %hopping to a quencher
                                    %Subtract exciton from current position
                                    np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5))) = np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5)))-1;
                                    quenched(aliveidx(cnt5)) = 1;
                                    alive(aliveidx(cnt5)) = 0;
                                    decaytime(aliveidx(cnt5)) = t;
                                end
                            end %If exciton neighbor = 0, cannot move forward, exciton stays put.
                        case 3 %Z axis
                            if np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5))+1) ~= 0
                                if np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5))+1) >= 1 %If a chromophore is present
                                    %Subtract exciton off of current position
                                    np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5))) = np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5)))-1;
                                    np.exposz(aliveidx(cnt5)) = np.exposz(aliveidx(cnt5))+1; %move exciton one step forward
                                    %Add exciton to forward position
                                    np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5))) = np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5)))+1;
                                elseif np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5))+1) == -1 %hopping to a quencher
                                    %Subtract exciton from current position
                                    np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5))) = np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5)))-1;
                                    quenched(aliveidx(cnt5)) = 1;
                                    alive(aliveidx(cnt5)) = 0;
                                    decaytime(aliveidx(cnt5)) = t;
                                end
                            end %If exciton neighbor = 0, cannot move forward, exciton stays put.
                    end
                elseif exhop(cnt5) == -1 %exciton hops backward, if possible.
                    switch i
                        case 1 %X axis
                            if np.latpts(np.exposx(aliveidx(cnt5))-1, np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5))) ~= 0
                                if np.latpts(np.exposx(aliveidx(cnt5))-1, np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5))) >= 1 %If a chromophore is present
                                    %Subtract exciton off of current position
                                    np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5))) = np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5)))-1;
                                    np.exposx(aliveidx(cnt5)) = np.exposx(aliveidx(cnt5))-1; %move exciton one step backward
                                    %Add exciton to forward position
                                    np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5))) = np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5)))+1;
                                elseif np.latpts(np.exposx(aliveidx(cnt5))-1, np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5))) == -1 %hopping to a quencher
                                    %Subtract exciton from current position
                                    np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5))) = np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5)))-1;
                                    quenched(aliveidx(cnt5)) = 1;
                                    alive(aliveidx(cnt5)) = 0;
                                    decaytime(aliveidx(cnt5)) = t;
                                end
                            end %If exciton neighbor = 0, cannot move forward, exciton stays put.
                        case 2 %Y axis
                            if np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5))-1, np.exposz(aliveidx(cnt5))) ~= 0
                                if np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5))-1, np.exposz(aliveidx(cnt5))) >= 1 %If a chromophore is present
                                    %Subtract exciton off of current position
                                    np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5))) = np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5)))-1;
                                    np.exposy(aliveidx(cnt5)) = np.exposy(aliveidx(cnt5))-1; %move exciton one step forward
                                    %Add exciton to forward position
                                    np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5))) = np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5)))+1;
                                elseif np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5))-1, np.exposz(aliveidx(cnt5))) == -1 %hopping to a quencher
                                    %Subtract exciton from current position
                                    np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5))) = np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5)))-1;
                                    quenched(aliveidx(cnt5)) = 1;
                                    alive(aliveidx(cnt5)) = 0;
                                    decaytime(aliveidx(cnt5)) = t;
                                end
                            end %If exciton neighbor = 0, cannot move forward, exciton stays put.
                        case 3 %Z axis
                            if np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5))-1) ~= 0
                                if np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5))-1) >= 1 %If a chromophore is present
                                    %Subtract exciton off of current position
                                    np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5))) = np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5)))-1;
                                    np.exposz(aliveidx(cnt5)) = np.exposz(aliveidx(cnt5))-1; %move exciton one step forward
                                    %Add exciton to forward position
                                    np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5))) = np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5)))+1;
                                elseif np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5))-1) == -1 %hopping to a quencher
                                    %Subtract exciton from current position
                                    np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5))) = np.latpts(np.exposx(aliveidx(cnt5)), np.exposy(aliveidx(cnt5)), np.exposz(aliveidx(cnt5)))-1;
                                    quenched(aliveidx(cnt5)) = 1;
                                    alive(aliveidx(cnt5)) = 0;
                                    decaytime(aliveidx(cnt5)) = t;
                                end
                            end %If exciton neighbor = 0, cannot move forward, exciton stays put.
                    end %If result of probhop calc yields rand > prob, exciton does not move. Proceed without hopping.
                end
            end
            %Finished hopping in 1D, need to calc decay and FRET for all excitons post 1D hop
            
            %Radiative/Non-radiative decay:
            justdecayed = alive .* (exp(-param.dt/param.tau0) < rand([1 param.nex]));
            alive = alive-justdecayed;
            decayed = decayed+justdecayed;
            decayidx = find(justdecayed == 1);
            decaytime(decayidx) = t;

            %FRET
            if param.quench >= 1
                aliveidx = find(alive);
                nlive = sum(alive);
                rates = zeros([1 nlive]);
                %Calculate quencher radial distances relative to exciton:
                rsq = zeros([1 nlive]);
                idx = 1:length(np.qposx);
                for idx = 1:length(np.qposx), %Calculate radial distance of all excitons to each single quencher
                    dx = (np.qposx(idx)-np.exposx(aliveidx)).*param.latspc;
                    dy = (np.qposy(idx)-np.exposy(aliveidx)).*param.latspc;
                    dz = (np.qposz(idx)-np.exposz(aliveidx)).*param.latspc;
                    rsq = dx.*dx + dy.*dy + dz.*dz;
                end
                %calculate FRET probability prior to hopping:
                FRET = param.r0^6/param.tau0;
                rates = rates + FRET./(rsq.*rsq.*rsq);
                justquenched = exp(-rates.*param.dt) < rand([1 nlive]);
                alive(aliveidx) = alive(aliveidx)-justquenched;
                quenched(aliveidx) = quenched(aliveidx)+justquenched;
                quenchidx = find(justquenched == 1);
                decaytime(quenchidx) = t;
            end

        end %Finished hopping in 3D, increase time step
    end
    results.qeff(navg) = (param.nex-sum(alive)-sum(decayed))/param.nex;
    results.sumpop = results.sumpop+results.pop;
end
%Output quenching efficiency:
qeff = mean(results.qeff);
fprintf(1,'Quenching efficiency (1.0=max): %6.4f\n',qeff);

%Output mean exciton lifetime:
pop = results.sumpop;
fdecay = pop/sum(pop);
halflife = sum(fdecay.*results.time);
fprintf(1,'Mean exciton lifetime %7.1f ps\n',halflife/log(2));

%plot population as a function of time:
plot(results.time, fdecay,'o')
xlabel('Time (ps)')
ylabel('Exciton Population')
title('Simulated Intensity Decay')

%Estimate beta:
lt = log(results.time);
lf = -log(pop/max(pop));
idx=find(fdecay==max(fdecay));
if length(idx)>1, % check for the case of 2 or more points=max
    idx=max(idx);
end
lt=lt((idx+1):end); % remove some early time points
lf=lf((idx+1):end);
llf = log(lf);
idx=find(isinf(llf));
if length(idx)>0,
    minidx=min(idx);
    lt=lt(1:(minidx-1));
    llf=llf(1:(minidx-1));
end
pp=polyfit(lt,llf,1);
fprintf(1,'estimated beta (approx) %7.4f\n',pp(1));
