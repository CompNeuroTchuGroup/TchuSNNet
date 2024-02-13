
def allocate_pos_deltat(
    timesteps, offset, steps, presynaptic, postsynaptic, deltaT_steps, noSpikes
):
    counter = 0
    counter2 = 0
    flag = False
    for i in timesteps:
        if i < offset or noSpikes < 1:
            continue
        elif counter % steps == 0 and not flag:
            presynaptic.append(i)
            counter2 = 0
            flag = True

        elif counter2 % deltaT_steps == 0 and flag:
            postsynaptic.append(i)
            flag = False
            noSpikes -= 1
        counter2 += 1
        counter += 1


def allocate_neg_deltat(
    timesteps, offset, steps, presynaptic, postsynaptic, deltaT_steps, noSpikes
):
    counter = 0
    counter2 = 0
    flag = False
    for i in timesteps:
        if i < offset or noSpikes < 1:
            continue
        elif counter % steps == 0 and not flag:
            postsynaptic.append(i)
            counter2 = 0
            flag = True

        elif counter2 % deltaT_steps == 0 and flag:
            presynaptic.append(i)
            flag = False
            noSpikes -= 1
        counter2 += 1
        counter += 1


def allocate_zero_dt(
    timesteps, offset, steps, presynaptic, postsynaptic, deltaT_steps, noSpikes
):
    counter = 0
    for i in timesteps:
        if i < offset or noSpikes < 1:
            continue
        elif counter % steps == 0:
            postsynaptic.append(i)
            presynaptic.append(i)
            noSpikes -= 1

        counter += 1
