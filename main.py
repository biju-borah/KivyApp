from turtle import width
import math
from kivy.app import App
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.textinput import TextInput
from kivy.uix.button import Button
from kivy.uix.screenmanager import ScreenManager, Screen
from kivy.uix.label import Label
import pandas as pd

from bdb import effective
from tkinter import FLAT

value = {'phi': 0.75, 'lamda': 1}


def beta1_value():
    if value['fc'] <= 28:
        beta1 = 0.85
    else:
        beta1 = max(0.85 - (0.05*((value['fc'] - 28)/7)), 0.65)
    return beta1


def steel_areas():
    import math
    As = (math.pi) * (int(value['dbar']) ** 2) / 4 * int(value['n'])
    Asprime = (math.pi) * (int(value['dbarprime'])
                           ** 2) / 4 * (int(value['nprime']))
    return As, Asprime


def effective_depth():
    num_layer = next(
        layer for layer in range(1, 6) if (value["b"] - 2 * value['cc'] - 2 * value['tie'] - (int(value['n']) / layer) * int(value['dbar'])) / ((int(value['n']) / layer) - 1) > int(value['dbar']))
    mult = num_layer - 0.5
    deff = value['h'] - value['cc'] - value['tie'] - int(value['dbar'])*mult
    return deff, num_layer


def topbar_depth():
    num_layer = next(layer for layer in range(1, 6) if
                     (value["b"] - 2 * value['cc'] - 2 * value['tie'] - (int(value['nprime']) / layer) * int(value['dbarprime'])) / ((int(value['nprime']) / layer) - 1) > int(value['dbarprime']))
    mult = num_layer - 0.5
    dprime = value['cc'] + value['tie'] + int(value['dbarprime']) * mult
    return dprime


def stress_strain():
    import math
    A = 0.85 * value['fc'] * value["b"] * value['beta1']
    B = value['Asprime'] * value['ec'] * \
        value['Es'] - value['As'] * value['fy']
    C = -value['Asprime'] * value['ec'] * value['Es'] * value['dprime']
    D = (B ** 2) - (4 * A * C)

    c = (-B + math.sqrt(D)) / (2 * A)
    a = value['beta1'] * c

    es = value['ec'] * ((value['deff'] - c) / c)
    fs = es * es

    esprime = value['ec'] * ((c - value['dprime']) / c)
    fsprime = esprime * value['Es']

    return fs, fsprime, es, a, c


def strength_factor_classification():
    phi = 0.65 + (0.25 * ((value['es'] - value['ey']) / (0.005 - value['ey'])))
    if phi >= 0.90:
        phi, classify = 0.90, 'Under-reinforced'
    elif phi <= 0.65:
        phi, classify = 0.65, 'Over-reinforced, kindly recheck the value and try again'
    else:
        phi, classify = phi, 'Balanced'
    return phi, classify


def forces_capacity():
    C = (0.85 * value['fc'] * value['a'] * value["b"]) / 1000
    Cprime = (value['Asprime'] * (min(value['fy'], value['fsprime'])) / 1000)
    T = (value['As'] * (min(value['fy'], value['fs']))) / 1000

    Mn = (C * (value['deff'] - (value['a'] / 2))) + \
        (Cprime * (value['deff'] - value['dprime']))
    Mu = (value['phi'] * Mn) / 1000

    return Mu, T, C, Cprime


def steel_ratio_check():
    rho = (value['As'] - value['Asprime']) / (value["b"] * value['deff'])
    rho_max = (3 / 7) * ((0.85 * value['fc'] * value['beta1']) / value['fy'])
    rho_min1 = 1.4 / value['fy']
    rho_min2 = (value['fc'] ** (1 / 2)) / (4 * value['fy'])

    if rho <= rho_max and rho >= min(rho_min1, rho_min2):
        return 'Steel ratio is WITHIN limits: ' '\u03C1 = {: .4f}'.format(rho)

    elif rho < min(rho_min1, rho_min2):
        return 'Steel ratio is LESS THAN the minimum limit: ''\u03A1' + 'min={:.4f}'.format(rho)

    else:
        return 'Steel ratio EXCEEDS the maximum limit, Overly-reinforced: ' '\u03A1' + 'max={:.4f}'.format(rho)


def rebar_suggestion():
    dbar_sizes = [16, 20, 25, 28, 32]
    n_list = []
    num_layer_list = []

    dbarprime_sizes = []
    nprime_list = []

    for dbar in dbar_sizes:
        import math
        Abar = (math.pi) * (dbar ** 2) / 4

        rho_max = (3 / 7) * (0.85 * value['fc'] * value['beta1']) / value['fy']
        deff_initial = value['h'] - value['cc'] - value['tie'] - dbar * 0.5
        As_max = rho_max * value['b'] * deff_initial

        n = math.ceil(As_max / Abar)
        n_list.append(n)

        num_layer = next(
            layer for layer in range(1, 6) if (value['b'] - 2 * value['cc'] - 2 * value['tie'] - (n / layer) * dbar) / ((n / layer) - 1) > dbar)
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

    df = pd.DataFrame(dict)
    return df


def cal1():
    try:
        # print("Please enter material and section properties:")
        # fc = float(input("Concrete strength, f'c(MPa) = "))
        # fy = float(input("Steel strength, fy (MPa) = "))
        # b = float(input("Beam width, b (mm) = "))
        # h = float(input("Beam total depth, h (mm) = "))
        fc = int(value["Concrete Strength"])
        fy = int(value["Steel Strength"])
        b = int(value["Beam Width"])
        h = int(value["Beam Depth"])
        value['fc'] = fc
        value['fy'] = fy
        value['b'] = b
        value['h'] = h
        print(
            "\nThis program assumes: Concrete cover = 40mm, Stirrups/Links Diameter = 10mm")
    except ValueError:
        print("DATA INPUT ERROR! Try again ")
    else:
        Es = 200000
        value['Es'] = Es
        value['ec'] = 0.003
        value['ey'] = fy / Es
        value['cc'] = 40
        value['tie'] = 10

        value['beta1'] = beta1_value()
        df = rebar_suggestion()
        capacity = []
        classify_member = []
        steel_ratio = []

        for i in range(0, len(df.index)):
            value['dbar'], value['n'], value['num_layer'], value['dbarprime'], value['nprime'] = df.iloc[i, :]

            value['As'], value['Asprime'] = steel_areas()
            value['deff'], value['num_layer'] = effective_depth()
            value['dprime'] = topbar_depth()
            value['fs'], value['fsprime'], value['es'], value['a'], value['c'] = stress_strain()
            value['phi'], classify = strength_factor_classification()
            Mu, value['T'], value['C'], value['Cprime'] = forces_capacity()

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

        # time.sleep(1)
        print(df[['Bottom Bars', 'Qty(bottom)', 'Classification', 'Steel Ratio']])

        return f"\n-----------------------------------\nSuggested Reinforcements:\n----------------------------------\n{df[['Bottom Bars', 'Qty(bottom)', 'Classification', 'Steel Ratio']]}"


def cal2():
    answer = 'y'
    while answer.lower().startswith("y"):
        try:
            dbar = int(value["dbar"])
            n = int(value["n"])
            dbarprime = int(value["dbarprime"])
            nprime = int(value["nprime"])
        except ValueError:

            print("DATA INPUT ERROR! Try again ")
        if dbar < 10 or n < 2 or dbarprime < 10 or nprime < 2:
            print("DATA INPUT ERROR! Minimum 2nos and Minimum 10mm dia")
        else:
            value['As'], value['Asprime'] = steel_areas()
            value['deff'], value['num_layer'] = effective_depth()
            value['dprime'] = topbar_depth()
            value['fs'], value['fsprime'], value['es'], value['a'], value['c'] = stress_strain()
            value['phi'], classify = strength_factor_classification()
            Mu, value['T'], value['C'], value['Cprime'] = forces_capacity()

            print(
                f'\nBeam Capacity: Mu = {round(Mu)} kN-m, \nBeam is under {classify}, \n{steel_ratio_check()}')
            print("\nHave a good day!\n")
            return f'\nBeam Capacity: Mu = {round(Mu )} kN-m, \nBeam is under {classify}, \n{steel_ratio_check()}'
            # exit()
    else:
        print("\nHave a good day!\n")
        return '\nHave a good day!\n'
        # exit()

################### app ##############################


class MyForm(BoxLayout):
    def __init__(self, **kwargs):
        super(MyForm, self).__init__(**kwargs)
        self.orientation = 'vertical'

        # Add input fields
        self.concrete_strength = TextInput(
            hint_text="Concrete strength, f'c(MPa)")
        self.add_widget(self.concrete_strength)

        self.steel_strength = TextInput(
            hint_text="Steel strength, fy (MPa)")
        self.add_widget(self.steel_strength)

        self.beam_width = TextInput(
            hint_text="Beam width, b (mm)")
        self.add_widget(self.beam_width)

        self.beam_depth = TextInput(
            hint_text="Beam total depth, h (mm)")
        self.add_widget(self.beam_depth)

        # Add submit button
        self.submit_button = Button(
            text='Submit', size_hint=(1.0, 1.0), size=(100, 50))
        self.submit_button.bind(on_press=self.submit_form)
        self.add_widget(self.submit_button)

    def submit_form(self, instance):
        # Code to handle form submission
        con_stren = self.concrete_strength.text
        steel_stren = self.steel_strength.text
        bw = self.beam_width.text
        bd = self.beam_depth.text
        value["Concrete Strength"] = con_stren
        value["Steel Strength"] = steel_stren
        value["Beam Width"] = bw
        value["Beam Depth"] = bd
        text = cal1()
        self.remove_widget(self.concrete_strength)
        self.remove_widget(self.steel_strength)
        self.remove_widget(self.beam_width)
        self.remove_widget(self.beam_depth)
        self.remove_widget(self.submit_button)
        # self.remove_widget(self.submit_button)
        self.label = Label(text=text)
        self.label2 = Label(
            text="\nThis program assumes: Concrete cover = 40mm, Stirrups/Links Diameter = 10mm", size_hint=(1.0, 0.2))
        self.continue_btn = Button(
            text='Continue', size_hint=(1.0, 0.5), size=(100, 50))
        self.continue_btn.bind(on_press=self.submit)
        self.add_widget(self.continue_btn)
        self.add_widget(self.label2)
        self.add_widget(self.label)

        # app = App.get_running_app()
        # app.root.current = 'second'
        # mainfun()

    def submit(self, instance):
        app = App.get_running_app()
        app.root.current = 'second'


class MyForm2(BoxLayout):
    def __init__(self, **kwargs):
        super(MyForm2, self).__init__(**kwargs)
        self.orientation = 'vertical'

        # Add input fields
        self.br = TextInput(
            hint_text="Bottom Rebar dia(mm)", height=20)
        self.add_widget(self.br)

        self.brt = TextInput(
            hint_text="Bottom Rebar Total qty (no.)", height=20)
        self.add_widget(self.brt)

        self.tr = TextInput(
            hint_text="Top Rebar dia(mm)", height=20)
        self.add_widget(self.tr)

        self.trt = TextInput(
            hint_text="Top Rebar Total qty (no.)", height=20)
        self.add_widget(self.trt)

        # Add submit button
        self.back_btn = Button(
            text='Back', size_hint=(1.0, 1.0), size=(100, 50))
        self.back_btn.bind(on_press=self.back2_btn)
        self.add_widget(self.back_btn)

        self.submit_button = Button(
            text='Submit', size_hint=(1.0, 1.0), size=(100, 50))
        self.submit_button.bind(on_press=self.submit_form)
        self.add_widget(self.submit_button)

    def back2_btn(self, instance):
        app = App.get_running_app()
        app.root.current = 'main'

    def back(self, instance):
        self.add_widget(self.br)
        self.add_widget(self.brt)
        self.add_widget(self.tr)
        self.add_widget(self.trt)
        self.add_widget(self.back_btn)
        self.add_widget(self.submit_button)
        self.remove_widget(self.label)
        self.remove_widget(self.formula)
        self.remove_widget(self.back_button)
        self.remove_widget(self.go_to_first_btn)

    def go_to_first(self, instance):
        screen2 = App.get_running_app().root.get_screen('second')
        screen2.reset()
        app = App.get_running_app()
        app.root.current = 'main'

    def submit_form(self, instance):
        # Code to handle form submission
        br = self.br.text
        brt = self.brt.text
        tr = self.tr.text
        trt = self.trt.text
        value["dbar"] = br
        value["n"] = brt
        value["dbarprime"] = tr
        value["nprime"] = trt
        text = cal2()
        print(text)
        self.remove_widget(self.br)
        self.remove_widget(self.brt)
        self.remove_widget(self.tr)
        self.remove_widget(self.trt)
        self.remove_widget(self.submit_button)
        self.remove_widget(self.back_btn)
        # φ * 0.17 * λ * √(f` c)*b w*d
        shear_strength = value['phi'] * 0.17 * value['lamda'] * \
            math.sqrt(int(value['fc'])) * value['b'] * value['h']
        self.formula = Label(
            text=f'Maximum steel ratio for beam = 0.04bD\nMinimum steal ratio for beam = 0.002bD')
        self.add_widget(self.formula)
        self.label = Label(
            text=f'Beam Shear Capacity : {int(shear_strength/1000)} KN/m2\n{text}')
        self.add_widget(self.label)
        self.back_button = Button(
            text='Back', size_hint=(1.0, 0.5), size=(100, 50))
        self.back_button.bind(on_press=self.back)
        self.add_widget(self.back_button)
        self.go_to_first_btn = Button(
            text='Reset', size_hint=(1.0, 0.5), size=(100, 50))
        self.go_to_first_btn.bind(on_press=self.go_to_first)
        self.add_widget(self.go_to_first_btn)
        screen1 = App.get_running_app().root.get_screen('main')
        screen1.reset()
        # app = App.get_running_app()
        # app.root.current = 'third'


class MainScreen(Screen):
    def __init__(self, **kwargs):
        super(MainScreen, self).__init__(**kwargs)
        # self.add_widget(Label(text="Welcome to the Main Screen"))
        # button = Button(text="Go to Second Screen")
        # button.bind(on_press=self.go_to_second_screen)
        # self.add_widget(button)
        self.form = MyForm()
        self.add_widget(self.form)

    # def reset_screen(self):
    #     self.reset_widget(self.form)

    def go_to_second_screen(self, instance):
        app = App.get_running_app()
        app.root.current = 'second'

    def reset(self):
        self.remove_widget(self.form)
        new_form = MyForm()
        self.add_widget(new_form)


class SecondScreen(Screen):
    def __init__(self, **kwargs):
        super(SecondScreen, self).__init__(**kwargs)
        self.form = MyForm2()
        self.add_widget(self.form)

    def go_to_main_screen(self, instance):
        app = App.get_running_app()
        app.root.current = 'main'

    def reset(self):
        self.remove_widget(self.form)
        new_form = MyForm2()
        self.add_widget(new_form)


class MyScreenManager(ScreenManager):
    pass


class MyApp(App):
    def build(self):
        # return MyForm()
        sm = MyScreenManager()
        sm.add_widget(MainScreen(name='main'))
        sm.add_widget(SecondScreen(name='second'))
        return sm


if __name__ == '__main__':
    MyApp().run()
