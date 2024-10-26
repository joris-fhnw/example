'''This file creates an object Wood and has a number of instances relevant for
wood combustion at 1 bar total pressure'''

import numpy as np
from scipy.optimize import minimize
from scipy import optimize
# import IG_Props_Comb_ed2 as Props
import CoolProp.CoolProp as CP
import cantera as ct
import matplotlib.pyplot as plt

# MW, Hf, C, H, O, N, S
# species of the substitute fuel : need to add N2 here...
H2 = np.array([2.01588, 0, 0, 2, 0, 0, 0])  # 0
H2O = np.array([18.0153, -241.825, 0, 2, 1, 0, 0])  # 1
CO = np.array([28.0101, -110.535, 1, 0, 1, 0, 0])  # 2
CO2 = np.array([44.01, -393.508, 1, 0, 2, 0, 0])  # 3
CH4 = np.array([16.042, -74.6, 1, 4, 0, 0, 0])  # 4
C2H6 = np.array([30.069, -83.851, 2, 6, 0, 0, 0])  # 5
C3H8 = np.array([44.096, -104.68, 3, 8, 0, 0, 0])  # 6
SO2 = np.array([64.064, -296.81, 0, 0, 2, 0, 1])  # 7
NO2 = np.array([46.0055, 33.10, 0, 0, 2, 1, 0])  # 8

zeta_O2_Air = 0.2314
zeta_N2_Air = 1 - zeta_O2_Air


class Wood(object):
    def __init__(self, fuel_type, gamma, water_content, Hu_dry=0):
        """
        Pass parameters describing molecules
        """
        # ! name
        self.type = fuel_type
        # extracting the masses of the elements in the dry wood
        self.C = 0
        self.H = 0
        self.O = 0
        self.N = 0
        self.S = 0
        for i in range(len(gamma)):
            if gamma[i] == 'C':
                self.C = gamma[i + 1]
            if gamma[i] == 'H':
                self.H = gamma[i + 1]
            if gamma[i] == 'O':
                self.O = gamma[i + 1]
            if gamma[i] == 'N':
                self.N = gamma[i + 1]
            if gamma[i] == 'S':
                self.S = gamma[i + 1]

        # normalizing these values to obtain mass fractions of dry wood
        ElTotal = self.C + self.H + self.O + self.N + self.S
        self.C = self.C / ElTotal
        self.H = self.H / ElTotal
        self.O = self.O / ElTotal
        self.N = self.N / ElTotal
        self.S = self.S / ElTotal

        self.w = water_content
        self.u = self.w / (1 - self.w)  # 1 kg of dry wood has u kg water
        self.Hu_dry = Hu_dry

        # normalized mass fractions of wet wood
        Total_w = self.C + self.H + self.O + self.N + self.S + self.u
        self.C_w = self.C / Total_w
        self.H_w = (self.H + self.u * (2 * 1.00794 / 18.0149)) / Total_w
        self.O_w = (self.O + self.u * (15.994 / 18.0149)) / Total_w
        self.N_w = self.N / Total_w
        self.S_w = self.S / Total_w
        self.Total_w = Total_w

        # molar composition of wet wood, basis 1 kg dry wood
        Cm_w = self.C / 12.0107
        Hm_w = (self.H + self.u * (2 * 1.00794 / 18.0149)) / 1.00794
        Om_w = (self.O + self.u * (15.994 / 18.0149)) / 15.994
        Nm_w = self.N / 14.0067
        Sm_w = self.S / 32.0
        ElTotalm = Cm_w + Hm_w + Om_w + Nm_w + Sm_w

        # Molar composition of wet wood, basis is C1;
        self.Cm_w_norm = Cm_w / ElTotalm
        self.Cm_w = 1.0
        self.Hm_w = (Hm_w / ElTotalm) / self.Cm_w_norm
        self.Om_w = (Om_w / ElTotalm) / self.Cm_w_norm
        self.Nm_w = (Nm_w / ElTotalm) / self.Cm_w_norm
        self.Sm_w = (Sm_w / ElTotalm) / self.Cm_w_norm
        self.mw = self.Cm_w * 12.0107 + self.Hm_w * 1.00794 + self.Om_w * 15.994 + \
                  self.Nm_w * 14.0067 + self.Sm_w * 32.0

        # Coefficients for balancing the molar reaction for complete oxidation of wet wood
        # C(Cm)H(Hm)S(Sm)O(Om)N(Nm) + lambda*a(O2 + 3.76N2) = bCO2 + cH2O + dSO2 + (lambda*e + f)N2 + (lambda-1)aO2
        self.a = (2 * self.Cm_w + self.Hm_w / 2 + 2 * self.Sm_w - self.Om_w) / 2
        self.b = self.Cm_w
        self.c = self.Hm_w / 2
        self.d = self.Sm_w
        self.e = self.a * 3.76
        self.f = self.Nm_w / 2

        # print('The fuel is C_%0.1f' % (self.Cm_w), 'H_%0.3f' % (self.Hm_w), 'O_%0.3f' % (self.Om_w), \
        #       'with "Wassergehalt" of %0.2f' % (self.w))

    def omin(self):  # this is calculated for 1 kg wet wood
        return 2.664 * self.C_w + 7.937 * (self.H / self.Total_w) + 0.998 * self.S_w - (self.O / self.Total_w)

    def Omin(self):  # this is calculated for 1 kmol wet wood
        return self.omin() * self.mw / 31.988

    def lmin(self):  # this is calculated for 1 kg wet wood
        """
        this is calculated for 1 kg wet wood
        :return: lmin [kg air/ kg wet wood]
        """
        return self.omin() / zeta_O2_Air

    def lmin_dry(self):  # this calculates lmin for 1 kg dry wood (as in Thermo script)
        """
        this is calculated for 1 kg dry wood
        :return: lmin [kg air/ kg dry wood]
        """
        return self.lmin() * (1 + self.u)

    def Lmin(self):  # this is calculated for 1 kmol wet wood
        """
        this is calculated for 1 kmol wet wood
        :return: lmin [kmol air/ kmol wet wood]
        """

        return self.lmin() * self.mw / 29.871

    def Hu_tr(self):  # in MJ/kg dry wood
        if self.Hu_dry == 0:
            return 34.8 * self.C + 93.9 * self.H + 6.3 * self.N - 10.8 * self.O
        else:
            return self.Hu_dry

    def Hu_nass(self):  # in MJ/kg wet wood
        """
        lower heating value of the wet wood in MJ/kg wet wood
        :return: lower heating value in MJ/kg
        """
        return self.Hu_tr() * (1 - self.w) - 2.443 * self.w

    def Hf_nass(self):  # in MJ/kg wet wood
        # for fuels with S and N
        Hf_nassm = self.Cm_w * CO2[1] + (self.Hm_w) / 2 * H2O[1] + (self.Sm_w) / 1 * SO2[1] + (self.Nm_w) / 1 * NO2[
            1] + self.Hu_nass() * self.mw
        return Hf_nassm / self.mw


    def T_ad(self, Tair_in, Rel_humidity, Lambda):  # Adiabatic flame temperature calculation (full oxidation)
        # Basis is 1 kg of wet fuel
        def H_fcn(T):
            return self.H_out(Tair_in, T, Rel_humidity, Lambda) - self.H_in(Tair_in, Rel_humidity, Lambda)

        T = 1200  # Initial Guess for Tad in K
        return optimize.fsolve(H_fcn, T)

    def ydot_wood_wet(self,x_o2_dry):
        """
        Berechnet wie viel feuchter Brennstoff für 1 kg trockene Luft notwendig ist
        :param x_o2_dry: [-] -> [Volumenbruch] Sauerstoffgehalt im Abgas
        :return: [kg Brennstoff feucht/ kg Luft fe]
        """
        m_wood_per_air = self.lmin()*self.Lamb_O2_tr(x_o2_dry)  # [kg Brennstoff/ kg Luft]
        return m_wood_per_air

    def ydot_wood_dry(self,x_o2_dry):
        """
        Berechnet wie viel trockener Brennstoff für 1 kg trockene Luft notwendig ist
        :param x_o2_dry: [-] -> [Volumenbruch] Sauerstoffgehalt im Abgas
        :return: [kg Brennstoff feucht/ kg Luft fe]
        """
        m_wood_per_air = self.lmin_dry()*self.Lamb_O2_tr(x_o2_dry)  # [kg Brennstoff/ kg Luft]
        return m_wood_per_air

    def SubstFuel(self):  # determining the mass fractions of Substitute Fuel for 1 kg wet wood
        Hf = self.Hf_nass()  # heat of formation of 1 kg wet wood
        C = self.C_w  # wt. of C in 1 kg wet wood
        H = self.H_w  # wt. of H in 1 kg wet wood
        O = self.O_w

        def objective(x):
            return abs(x[0] * H2[1] / H2[0] + x[1] * H2O[1] / H2O[0] + x[2] * CO[1] / CO[0] + \
                       x[3] * CO2[1] / CO2[0] + x[4] * CH4[1] / CH4[0] + x[5] * C2H6[1] / C2H6[0] + \
                       x[6] * C3H8[1] / C3H8[0] - Hf)

        def constraint1(x):  # C Balance
            return 12.0107 * (x[2] * CO[2] / CO[0] + x[3] * CO2[2] / CO2[0] + x[4] * CH4[2] / CH4[0] + \
                              x[5] * C2H6[2] / C2H6[0] + x[6] * C3H8[2] / C3H8[0]) - C

        def constraint2(x):  # H Balance
            return 1.00794 * (x[0] * H2[3] / H2[0] + x[1] * H2O[3] / H2O[0] + \
                              x[4] * CH4[3] / CH4[0] + x[5] * C2H6[3] / C2H6[0] + x[6] * C3H8[3] / C3H8[0]) - H

        def constraint3(x):  # O Balance
            return 15.994 * (x[1] * H2O[4] / H2O[0] + x[2] * CO[4] / CO[0] + \
                             x[3] * CO2[4] / CO2[0]) - O

        # initial guesses
        n = 7
        x0 = np.zeros(7)
        x0[0] = 0.0
        x0[1] = 0.16
        x0[2] = 0.0
        x0[3] = 0.55
        x0[4] = 0.25
        x0[5] = 0.0
        x0[6] = 0.0

        # show initial objective
        # print('Initial Objective: ' + str(objective(x0)))

        # optimize
        b = (0.0, 1.0)
        bnds = (b, b, b, b, b, b, b)
        con1 = {'type': 'eq', 'fun': constraint1}
        con2 = {'type': 'eq', 'fun': constraint2}
        con3 = {'type': 'eq', 'fun': constraint3}

        cons = ([con1, con2, con3])
        solution = minimize(objective, x0, method='SLSQP', \
                            bounds=bnds, constraints=cons)
        x = solution.x
        # show final objective
        # print('Final Objective: ' + str(objective(x)))
        # print(sum(x))
        # print('con1',constraint1(x))
        # print('con2',constraint2(x))
        # print('con3',constraint3(x))
        if abs(constraint1(x)) > 1e-3:
            print('warning con1')
        if abs(constraint2(x)) > 1e-3:
            print('warning con2')
        if abs(constraint3(x)) > 1e-3:
            print('warning con3')
        if abs(objective(x)) > 1e-3:
            print('warning Hf')
        return x  # mass fractions for substitute fuel

    def SubstFuel_mol(self):  # determining the mole fractions of Substitute Fuel for wet wood
        massfcn = self.SubstFuel()
        TotalMol = massfcn[0] / H2[0] + massfcn[1] / H2O[0] + massfcn[2] / CO[0] + massfcn[3] / CO2[0] + \
                   massfcn[4] / CH4[0] + massfcn[5] / C2H6[0] + massfcn[6] / C3H8[0]

        y = np.zeros(7)
        y[0] = (massfcn[0] / H2[0]) / TotalMol
        y[1] = (massfcn[1] / H2O[0]) / TotalMol
        y[2] = (massfcn[2] / CO[0]) / TotalMol
        y[3] = (massfcn[3] / CO2[0]) / TotalMol
        y[4] = (massfcn[4] / CH4[0]) / TotalMol
        y[5] = (massfcn[5] / C2H6[0]) / TotalMol
        y[6] = (massfcn[6] / C3H8[0]) / TotalMol
        # print(sum(y))
        y_d = {'H2': y[0], 'H2O': y[1], 'CO': y[2], 'CO2': y[3], \
         'CH4': y[4], 'C2H6': y[5], 'C3H8': y[6]}

        return y,y_d  # moles fractions for substitute fuel

    def mw_SubstFuel(self):  # returns the molecular wt of substitute fuel
        y = self.SubstFuel_mol()
        return y[0] * H2[0] + y[1] * H2O[0] + y[2] * CO[0] + y[3] * CO2[0] + y[4] * CH4[0] + \
               y[5] * C2H6[0] + y[6] * C3H8[0]

    def T_ad_ct(self, Tair_in, Rel_humidity, Lambda):
        gas = ct.Solution('gri30.xml')
        nsp = gas.n_species

        p = 1  # 1 bar combustion pressure
        P_H2O_Air = Rel_humidity * CP.PropsSI("P", "T", Tair_in, "Q", 1, "water") / 10 ** 5  # in bar

        # find the needed species indices
        ih2 = gas.species_index('H2')
        ih2o = gas.species_index('H2O')
        ico = gas.species_index('CO')
        ico2 = gas.species_index('CO2')
        ich4 = gas.species_index('CH4')
        ic2h6 = gas.species_index('C2H6')
        ic3h8 = gas.species_index('C3H8')
        io2 = gas.species_index('O2')
        in2 = gas.species_index('N2')

        y = self.SubstFuel_mol()
        phi = 1 / Lambda
        # molar amounts of fuel and air
        x = np.zeros([nsp])
        x[io2] = y[0] / 2 + y[2] / 2 + y[4] * (2) + y[5] * (3.5) + y[6] * (
            5)  # this is Omin for the substitute fuel (mols O2 for 1 mol substitute fuel)
        x[in2] = 3.76 * x[io2]
        x_H2O_Air = P_H2O_Air / p * (x[io2] + x[in2]) / (1 - P_H2O_Air / p)

        x[ih2] = phi * y[0]
        x[ih2o] = phi * y[1] + x_H2O_Air
        x[ico] = phi * y[2]
        x[ico2] = phi * y[3]
        x[ich4] = phi * y[4]
        x[ic2h6] = phi * y[5]
        x[ic3h8] = phi * y[6]

        gas.TPX = Tair_in, 1e5, x
        gas.equilibrate('HP')
        return gas.T

    # Exhaust Gas composition (only for Lambda > 1), using the molar comp of wet wood
    def CO2(self, Lambda):
        return self.b / (self.b + self.c + self.d + (Lambda * self.e + self.f) + (Lambda - 1) * self.a)

    def CO2_tr(self, Lambda):
        return self.b / (self.b + self.d + (Lambda * self.e + self.f) + (Lambda - 1) * self.a)

    def H2O(self, Lambda):
        return self.c / (self.b + self.c + self.d + (Lambda * self.e + self.f) + (Lambda - 1) * self.a)

    def SO2(self, Lambda):
        return self.d / (self.b + self.c + self.d + (Lambda * self.e + self.f) + (Lambda - 1) * self.a)

    def N2(self, Lambda):
        return ((Lambda * self.e + self.f) * self.a) / (
                    self.b + self.c + self.d + (Lambda * self.e + self.f) + (Lambda - 1) * self.a)

    def O2(self, Lambda):
        return ((Lambda - 1) * self.a) / (self.b + self.c + self.d + (Lambda * self.e + self.f) + (Lambda - 1) * self.a)

    def O2_tr(self, Lambda):
        return ((Lambda - 1) * self.a) / (self.b + self.d + (Lambda * self.e + self.f) + (Lambda - 1) * self.a)

    def Lamb_CO2_tr(self, CO2_tr):
        def Lamb_fcn(Lamb):
            return self.CO2_tr(Lamb) - CO2_tr

        Lamb = 1.5  # intitial guess
        return optimize.fsolve(Lamb_fcn, Lamb)

    def Lamb_O2_tr(self, O2_tr):
        def Lamb_fcn(Lamb):
            return self.O2_tr(Lamb) - O2_tr

        Lamb = 1.5  # intitial guess
        return optimize.fsolve(Lamb_fcn, Lamb)

    def ctComb_eq_state(self, Tair_in, Rel_humidity, Lambda):
        eq_state = ct.Solution('gri30.xml')
        nsp = eq_state.n_species
        x = np.zeros([nsp, 1])

        p = 1  # 1 bar combustion pressure
        P_H2O_Air = Rel_humidity * CP.PropsSI("P", "T", Tair_in, "Q", 1, "water") / 10 ** 5  # in bar

        # find the needed species indices
        ih2 = eq_state.species_index('H2')
        ih2o = eq_state.species_index('H2O')
        ico = eq_state.species_index('CO')
        ico2 = eq_state.species_index('CO2')
        ich4 = eq_state.species_index('CH4')
        ic2h6 = eq_state.species_index('C2H6')
        ic3h8 = eq_state.species_index('C3H8')
        io2 = eq_state.species_index('O2')
        in2 = eq_state.species_index('N2')

        y = self.SubstFuel_mol()
        phi = 1 / Lambda
        # molar amounts of fuel and air
        x[io2, 0] = y[0] / 2 + y[2] / 2 + y[4] * (2) + y[5] * (3.5) + y[6] * (
            5)  # this is Omin for the substitute fuel (mols O2 for 1 mol substitute fuel)
        x[in2, 0] = 3.76 * x[io2, 0]
        x_H2O_Air = P_H2O_Air / p * (x[io2] + x[in2]) / (1 - P_H2O_Air / p)

        x[ih2, 0] = phi * y[0]
        x[ih2o, 0] = phi * y[1] + x_H2O_Air
        x[ico, 0] = phi * y[2]
        x[ico2, 0] = phi * y[3]
        x[ich4, 0] = phi * y[4]
        x[ic2h6, 0] = phi * y[5]
        x[ic3h8, 0] = phi * y[6]

        eq_state.TPX = Tair_in, 1e5, x
        eq_state.equilibrate('HP')
        print(eq_state.report())
        return eq_state

    def ctComb_eq(self, Tair_in, Rel_humidity):
        gas = ct.Solution('gri30.xml')
        nsp = gas.n_species

        p = 1  # 1 bar combustion pressure
        P_H2O_Air = Rel_humidity * CP.PropsSI("P", "T", Tair_in, "Q", 1, "water") / 10 ** 5  # in bar

        phi = np.zeros(100)
        Lambda = np.zeros(100)
        tad = np.zeros(100)
        xeq = np.zeros([nsp, 100])

        # find the needed species indices
        ih2 = gas.species_index('H2')
        ih2o = gas.species_index('H2O')
        ico = gas.species_index('CO')
        ico2 = gas.species_index('CO2')
        ich4 = gas.species_index('CH4')
        ic2h6 = gas.species_index('C2H6')
        ic3h8 = gas.species_index('C3H8')
        io2 = gas.species_index('O2')
        in2 = gas.species_index('N2')

        y = self.SubstFuel_mol()

        for i in range(0, 100):
            phi[i] = 1 + 0.09 * i
            Lambda[i] = 1 / phi[i]
            x = np.zeros([nsp, 1])
            # molar amounts of fuel and air
            x[io2, 0] = y[0] / 2 + y[2] / 2 + y[4] * (2) + y[5] * (3.5) + y[6] * (
                5)  # this is Omin for the substitute fuel (mols O2 for 1 mol substitute fuel)
            x[in2, 0] = 3.76 * x[io2, 0]
            x_H2O_Air = P_H2O_Air / p * (x[io2] + x[in2]) / (1 - P_H2O_Air / p)

            x[ih2, 0] = phi[i] * y[0]
            x[ih2o, 0] = phi[i] * y[1] + x_H2O_Air
            x[ico, 0] = phi[i] * y[2]
            x[ico2, 0] = phi[i] * y[3]
            x[ich4, 0] = phi[i] * y[4]
            x[ic2h6, 0] = phi[i] * y[5]
            x[ic3h8, 0] = phi[i] * y[6]

            gas.TPX = Tair_in, 1e5, x
            gas.equilibrate('HP')
            tad[i] = gas.T
            xeq[:, i] = gas.X

        total = plt.figure()
        diag1 = total.add_subplot(1, 3, 1)
        diag1.plot(Lambda, tad, label='Tad as Function of Lambda')
        diag1.set_xlabel('Lambda')
        diag1.set_ylabel('Adiabatic Flame Temperature, K')
        diag1.legend(loc='best')
        diag1.grid('on')

        diag2 = total.add_subplot(1, 3, 2)
        diag2.plot(Lambda, xeq[0, :], label='yH2')
        diag2.plot(Lambda, xeq[5, :], label='yH2O')
        diag2.plot(Lambda, xeq[13, :], label='yCH4')
        diag2.plot(Lambda, xeq[14, :], label='yCO')
        diag2.plot(Lambda, xeq[15, :], label='yCO2')
        diag2.plot(Lambda, xeq[47, :], label='yN2')
        diag2.set_xlabel('Lambda')
        diag2.set_ylabel('Molanteile')
        diag2.legend(loc='best')
        diag2.grid('on')

        diag3 = total.add_subplot(1, 3, 3)
        diag3.plot(tad, xeq[0, :], label='yH2')
        diag3.plot(tad, xeq[5, :], label='yH2O')
        diag3.plot(tad, xeq[13, :], label='yCH4')
        diag3.plot(tad, xeq[14, :], label='yCO')
        diag3.plot(tad, xeq[15, :], label='yCO2')
        diag3.plot(tad, xeq[47, :], label='yN2')
        diag3.set_xlabel('Adiabatic Flame Temperature, K')
        diag3.set_ylabel('Molanteile')
        diag3.legend(loc='best')
        diag3.grid('on')

        total.show()