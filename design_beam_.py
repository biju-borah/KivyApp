from bdb import effective
from tkinter import FLAT


def beta1_value():
    if fc <= 28:
        beta1 = 0.85
    else:
        beta1 = max (0.85 - (0.05*((fc - 28)/7)) , 0.65 )
    return beta1


def steel_areas():
    import math
    As = (math.pi) * (dbar ** 2) / 4 * n
    Asprime = (math.pi) * (dbarprime ** 2) / 4 * (nprime)
    return As, Asprime


def effective_depth():
    num_layer = next(
        layer for layer in range(1, 6) if (b - 2 * cc - 2 * tie - (n / layer) * dbar) / ((n / layer) - 1) > dbar)
    mult = num_layer - 0.5
    deff = h - cc - tie - dbar*mult
    return deff, num_layer


def topbar_depth():
    num_layer = next(layer for layer in range(1, 6) if
                     (b - 2 * cc - 2 * tie - (nprime / layer) * dbarprime) / ((nprime / layer) - 1) > dbarprime)
    mult = num_layer - 0.5
    dprime = cc + tie + dbarprime * mult
    return dprime


def stress_strain():
    import math
    A = 0.85 * fc * b * beta1
    B = Asprime * ec * Es - As * fy
    C = -Asprime * ec * Es * dprime
    D = (B ** 2) - (4 * A * C)

    c = (-B + math.sqrt(D)) / (2 * A)
    a = beta1 * c

    es = ec * ((deff - c) / c)
    fs = es * es

    esprime = ec * ((c - dprime) / c)
    fsprime = esprime * Es

    return fs, fsprime, es, a, c


def strength_factor_classification():
    phi = 0.65 + (0.25 * ((es - ey) / (0.005 - ey)))
    if phi >= 0.90:
        phi, classify = 0.90, 'Under-reinforced'
    elif phi <= 0.65:
        phi, classify = 0.65, 'Over-reinforced'
    else:
        phi, classify = phi, 'Balanced'
    return phi, classify


def forces_capacity():
    C = (0.85 * fc * a * b) / 1000
    Cprime = (Asprime * (min(fy, fsprime)) / 1000)
    T = (As * (min(fy, fs))) / 1000

    Mn = (C * (deff - (a / 2))) + (Cprime * (deff - dprime))
    Mu = (phi * Mn) / 1000

    return Mu, T, C, Cprime


def steel_ratio_check():
    rho = (As - Asprime) / (b * deff)
    rho_max = (3 / 7) * ((0.85 * fc * beta1) / fy)
    rho_min1 = 1.4 / fy
    rho_min2 = (fc ** (1 / 2)) / (4 * fy)

    if rho <= rho_max and rho >= min(rho_min1, rho_min2):
        return 'Steel ratio is WITHIN limits: ' '\u03C1 = {: .4f}'.format(rho)

    elif rho < min(rho_min1, rho_min2):
        return 'Steel ratio is LESS THAN the minimum limit: ''\u03A1' + 'min={:.4f}'.format(rho)

    else:
        return 'Steel ratio EXCEEDS the maximum limit,Overly-reinforced: ' '\u03A1' + 'max={:.4f}'.format(rho)


def rebar_suggestion():
    dbar_sizes = [20, 25, 28, 32]
    n_list = []
    num_layer_list = []

    dbarprime_sizes = []
    nprime_list = []

    for dbar in dbar_sizes:
        import math
        Abar = (math.pi) * (dbar ** 2) / 4

        rho_max = (3 / 7) * (0.85 * fc * beta1) / fy
        deff_initial = h - cc - tie - dbar * 0.5
        As_max = rho_max * b * deff_initial

        n = math.ceil(As_max / Abar)
        n_list.append(n)

        num_layer = next(
            layer for layer in range(1, 6) if (b - 2 * cc - 2 * tie - (n / layer) * dbar) / ((n / layer) - 1) > dbar)
        num_layer_list.append(num_layer)

        dbarprime = dbar
        dbarprime_sizes.append(dbarprime)

        if num_layer > 1 and n / num_layer > 2:
            nprime = 3
        else:
            nprime = 2
        nprime_list.append(nprime)

    dict = {'Bottom Bars': dbar_sizes, 'Qty(bottom)': n_list, 'Layer': num_layer_list, 'Top Bars': dbarprime_sizes,
            'Qty(top)': nprime_list}

    import pandas as pd
    df = pd.DataFrame(dict)
    return df


print('\n-------------------------------------------')
print("BEAM DESIGN CHECK - Rectangular Beam")
print('--------------------------------------------')

while True:
    try:
        print("Please enter material and section properties:")
        fc = float(input("Concrete strength, f'c(MPa) = "))
        fy = float(input("Steel strength, fy (MPa) = "))
        b = float(input("Beam width, b (mm) = "))
        h = float(input("Beam total depth, h (mm) = "))
        print("\nThis program assumes: Concrete cover = 40mm, Stirrups/Links Diameter = 10mm")
    except ValueError:
        print("DATA INPUT ERROR! Try again ")
    else:
        Es = 200000
        ec = 0.003
        ey = fy / Es
        cc = 40
        tie = 10

        beta1 = beta1_value()
        df = rebar_suggestion()
        capacity = []
        classify_member = []
        steel_ratio = []

        for i in range(0, len(df.index)):
            dbar, n, num_layer, dbarprime, nprime = df.iloc[i, :]

            As, Asprime = steel_areas()
            deff, num_layer = effective_depth()
            dprime = topbar_depth()
            fs, fsprime, es, a, c = stress_strain()
            phi, classify = strength_factor_classification()
            Mu, T, C, Cprime = forces_capacity()

            capacity.append(round(Mu))
            classify_member.append(classify)
            steel_ratio.append(steel_ratio_check())

        df['Capacity, Mu(kN-m)'] = capacity
        df['Classification'] = classify_member
        df['Steel Ratio'] = steel_ratio

        print('\n-----------------------------------')
        print('Suggested Reinforcements: ')
        print('----------------------------------')
        import time

        time.sleep(1)
        print(df)

        answer = input("\nWould you like to provide reinforcements to check? (y/n) ")
        while answer.lower().startswith("y"):
            try:
                dbar = int(input("Bottom Rebar dia(mm) = "))
                n = int(input("Bottom Rebar Total qty (no.) = "))
                dbarprime = int(input("Top Rebar dia(mm) = "))
                nprime = int(input("Top Rebar Total qty (no.) = "))
            except ValueError:

                print("DATA INPUT ERROR! Try again ")
            if dbar < 10 or n < 2 or dbarprime < 10 or nprime < 2:
                print("DATA INPUT ERROR! Minimum 2nos and Minimum 10mm dia")
            else:
                As, Asprime = steel_areas()
                deff, num_layer = effective_depth()
                dprime = topbar_depth()
                fs, fsprime, es, a, c = stress_strain()
                Mu, T, C, Cprime = forces_capacity()
                phi, classify = strength_factor_classification()

                print(f'\nBeam Capacity: Mu = {round(Mu )} kN-m, \nBeam is under {classify}, \n{steel_ratio_check()}')
                print("\nHave a good day!\n")
                exit()
        else:
            print("\nHave a good day!\n")
            exit()
