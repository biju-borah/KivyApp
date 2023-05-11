from kivy.app import App
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.textinput import TextInput
from kivy.uix.button import Button


from bdb import effective
from tkinter import FLAT


def beta1_value():
    if fc <= 28:
        beta1 = 0.85
    else:
        beta1 = max(0.85 - (0.05*((fc - 28)/7)), 0.65)
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
            text='Submit', size_hint=(None, None), size=(100, 50))
        self.submit_button.bind(on_press=self.submit_form)
        self.add_widget(self.submit_button)

    def submit_form(self, instance):
        # Code to handle form submission
        con_stren = self.concrete_strength.text
        steel_stren = self.steel_strength.text
        bw = self.beam_width.text
        bd = self.beam_depth.text
        print(
            f"Concrete Strength: {con_stren}, Steel Strength: {steel_stren}, Beam width: {bw}, Beam Depth: {bd}")


class MyApp(App):
    def build(self):
        return MyForm()


if __name__ == '__main__':
    MyApp().run()
