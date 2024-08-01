import sympy as sp
import numpy as np
import itertools
import matplotlib.pyplot as plt

class TransferFunction:
    '''
    This class is used when using the get_symbolic_transfert_function() method of the Circuit class,
    in order to extract the symbolic transfer function of a circuit object, between two nodes. 
    The object created must contain all the components attributes of the given circuit associated 
    (in order to do the substitutions of the components values in the symbolic transfer function).

    Normally, you should not need to use this class directly, but only through the 
    get_symbolic_transfert_function() method of the Circuit class, that will return an object of 
    this class.
    '''

    def __init__(self, H, components):
        self.sympyExpr = H
        self.components = components
        self.b, self.a = None, None  #Symbolic analog filter coefficients
        self.Nb, self.Na = None, None # Length of coefficients

    def sym_coeffs(self, norm=False):
        '''
        Returns the transfer function's symbolic coefficients.
        
        Parameters
        ----------
        norm : {True, False}, optional
            Normalize the coefficients w.r.t a[0].
        '''
        if self.b is None or self.a is None:
            num, denum = sp.fraction(self.sympyExpr)
            self.b = sp.Poly(num, sp.symbols('s')).all_coeffs()
            self.a = sp.Poly(denum, sp.symbols('s')).all_coeffs()
            
        # Normalize coeffs
        if norm == True:
            self.__normalized()
        
        # Update coeffs size
        self.Nb = len(self.b)
        self.Na = len(self.a)
            
        return self.b, self.a

    def __normalized(self):
        '''
        Normalize the transfer function to have self.a[0] = 1.
        '''

        a0 = self.a[0]
        self.a = [coeff/a0 for coeff in self.a]
        self.b = [coeff/a0 for coeff in self.b]

        self.sympyExpr = sp.Poly(self.b, sp.symbols('s')) / sp.Poly(self.a, sp.symbols('s'))

    def num_coeffs(self, component_values=None, combination='nested'):
        """
        Return the numerical coefficients `b_num` and `a_num` of the analog filter transfer function.
        The coefficients are calculated by substituting the component values in the symbolic transfer function.

        Parameters
        ----------
        component_values : {None, dict}, optional
            A dictionary of component values. The keys are the component symbols, and the values are the component values.
            * If component_values is None, the default values of the components set in the Circuit object will be used.
            
            * If component_values is a dictionary with one key and a 1D array key values, 
            the function will return the a and b numerical coefficients for each value of the array, in a 2D array.
            
            * If component_values is a dictionary with multiple keys and with their associated 1D array key values,
            the function will return the a and b numerical coefficients for each combination of values in a (N+1)-D array,
            where N is the number of keys in the dictionary.

        combination : {'nested', 'parallel'}, optional
            * If the `component_values` dictionary has multiple keys, the 'combinations' argument specifies how to combine the values.
            * If the `component_values` dictionary has only one key, the 'combinations' argument is ignored 
            (for the moment not really, but it should works with 'nested' or 'parallel' for one key too, even if it's not optimal)

            * If `'nested'`, return all the possible combinations of the component values. 
            For example, if the dictionary is `{'R1': [1, 2], 'R2': [3, 4, 5]}`, 
            the `'nested'` option will return the combinations `[[(1, 3), (1, 4), (1, 5)], [(2, 3), (2, 4), (2,5)]]`.
            * If 'parallel', return only the combinations of the component values in the order of the dictionary keys.
            Therefore, the same number of values must be provided for each key in the dictionary.
            For example, if the dictionary is `{'R1': [1, 2], 'R2': [3, 4]}`, 
            the 'parallel' option will return the combinations `[(1, 3), (2, 4)]`.

        Returns
        -------
        b_num : list
            The numerator coefficients of the transfer function.
        a_num : list
            The denominator coefficients of the transfer function.
        """

        if self.b is None or self.a is None:
            self.b, self.a = self.sym_coeffs(norm=True)
        
        if component_values is None:
            component_values = {component.symbol: component.value for component in self.components}

            b_num = np.array([float(coeff.subs(component_values)) for coeff in self.b])
            a_num = np.array([float(coeff.subs(component_values)) for coeff in self.a])

            return b_num, a_num

        keys, values = zip(*component_values.items())

        if combination == 'nested':
            combinations = list(itertools.product(*values))
            shape = [len(values_set) for values_set in values] + [len(self.b)]
        elif combination == 'parallel':
            if not all(len(value) == len(values[0]) for value in values):
                raise ValueError("All values in the component_values dictionary must have the same length if 'parallel' is selected.")
            combinations = list(zip(*values))
            shape = (len(combinations), len(self.b))
        else:
            raise ValueError("The 'combinations' argument must be either 'nested' or 'parallel'.")

        b_num, a_num = [], []
        for comb in combinations:
            substitutions = { component.symbol: comb[keys.index(str(component.symbol))]
                            if str(component.symbol) in keys else component.value for component in self.components}

            b_num_temp = [float(coeff.subs(substitutions)) for coeff in self.b]
            a_num_temp = [float(coeff.subs(substitutions)) for coeff in self.a]

            b_num.append(b_num_temp)
            a_num.append(a_num_temp)

        return np.array(b_num).reshape(shape), np.array(a_num).reshape(shape)


    def coeffs(self, values = 'symb', scheme = None, Fs = None, norm = False, combination = 'nested', simple=False):
        """
        Return b & a coefficients of the transfer function.

        Parameters
        ----------
        - values : 'symb', 'num', or dict like {'R1': [1e3, 2e6], 'R2': [3e2, 4e2]}
            * If  'sym', return the symbolic coefficients.
            * If  'num', return the numerical coefficients for the default values of the components set in the Circuit object.
            * If a dictionary is provided, return the numerical coefficients for the given values.
        - scheme : {'frwrd', 'bckwrd', 'blnr'}, or None, optional
            The desired discretization scheme.
            * If None, the function will return the analog filter coefficients.
        - Fs : float, optional
            The sampling frequency.
        - norm (bool): Returns normalized coefficients w.r.t a[0]
        - combination : {'nested', 'parallel'}, optional
            * If the `component_values` dictionary has multiple keys, the 'combinations' argument specifies how to combine the values.
            * If the `component_values` dictionary has only one key, the 'combinations' argument is ignored
        - simple (bool): Analog coefficients are expressed implicitly, returns generic discretized coefficients.
        Returns
        -------
        b, a : np.array
            The b & a coefficients of the transfer function.
        """
        
        # First, we derive the analog filter coefficients:

        if values == 'symb':
            b, a = self.sym_coeffs()
        elif values == 'num' or isinstance(values, dict):
            component_values = None if values == 'num' else values
            b, a = self.num_coeffs(component_values=component_values, combination=combination)
        else:
            raise ValueError("The 'values' argument must be either 'symb', 'num', or a dictionary of component values.")

        
        # Then we use the eval() function to derive the digital filter coefficients (if the user want) 
        # corresponding to the analog filter coefficients, either they are symbolic or numerical.

        if scheme is not None:
            b, a = self.coeffz(scheme=scheme, Fs=Fs, norm=norm, simple=simple)
        return b, a

    def sym_coeffz(self, scheme='blnr', norm=False, simple=False):
        """
        Returns the symbolic discretized coefficients for a given order N & discretization scheme.

        Parameters:
        - scheme (str): Desired discretization scheme ('frwrd', 'bckwrd', 'blnr').
        - norm (bool): Returns normalized coefficients w.r.t a[0].
        - simple (bool): Analog coefficients are expressed implicitly, returns generic discretized coefficients.

        Returns:
        - Bd (np.array): Discretized numerator coefficients.
        - Ad (np.array): Discretized denominator coefficients.
        """
        s, z, Ts = sp.symbols('s z T_s')

        #create a dictionary of the different discretization schemes
        #therefore, additonal schemes can be added easily
        schemes = {
            'frwrd': (z - 1) / Ts,
            'bckwrd': (z - 1) / (Ts * z),
            'blnr': 2 / Ts * (z - 1) / (z + 1)
        }

        if scheme not in schemes:
            raise ValueError("Invalid scheme given")
        
        # Make sure coeffs exist
        b, a = self.sym_coeffs()
        
        if simple == True:
            # Generalize analog coeffs
            b = sp.symbols(f'b_0:{self.Nb}')[::-1]
            a = sp.symbols(f'a_0:{self.Na}')[::-1]
        
        # Set discretization scheme
        ds = schemes[scheme]
        
        # Create variables
        B = sp.Poly(b, s)
        A = sp.Poly(a, s)

        # Replace s for chosen method
        B = B.as_expr().subs(s, ds) * (z + 1)**(self.Nb - 1)
        A = A.as_expr().subs(s, ds) * (z + 1)**(self.Na - 1)

        # Simplify
        H = B / A
        B, A = sp.fraction(H.simplify())

        # Group terms
        B = sp.collect(sp.expand(B / z**(self.Nb - 1) ), 1 / z)
        A = sp.collect(sp.expand(A / z**(self.Na - 1) ), 1 / z)

        # Store coefficients in order in a list
        Bd = np.zeros(self.Nb, dtype=object)
        Ad = np.zeros(self.Na, dtype=object)
        
        # Extract coeffs
        for n in range(self.Nb):
            Bd[n] = B.coeff(z, -n)
            
        for n in range(self.Na):
            Ad[n] = A.coeff(z, -n)

        return Bd, Ad



    def num_coeffz(self, scheme='blnr', Fs=None):
        """
        Returns the symbolic or real discretized coefficients for a given array of analog coefficients.

        Parameters:
        - scheme (str): Discretization scheme ('frwrd', 'bckwrd', 'blnr').
        - Fs (float): Samplerate.

        Returns:
        - Bd (np.array): Discretized numerator coefficients.
        - Ad (np.array): Discretized denominator coefficients.
        """
        
        B, A = self.sym_coeffz(scheme, simple=True)
        
        if Fs == None:
            raise ValueError('Fs was not given.')
        subs_dict = {'T_s': 1. / Fs}

        b, a = self.num_coeffs()
        
        # Reverse coeffs
        b = b[::-1]
        a = a[::-1]
        
        for n in range(self.Nb):
            subs_dict[f'b_{n}'] = b[n]
        
        for n in range(self.Na):
            subs_dict[f'a_{n}'] = a[n]

        Bd = np.array([B[n].subs(subs_dict) for n in range(self.Nb)])
        Ad = np.array([A[n].subs(subs_dict) for n in range(self.Na)])
        
        Bd = Bd / Ad[0]
        Ad = Ad / Ad[0]

        return Bd.astype(np.float64), Ad.astype(np.float64)



    def coeffz(self, values='sym', scheme='blnr', Fs=None, norm=False, simple=False):
        """
        Returns the symbolic or real discretized coefficients for a given array of analog coefficients.

        Parameters:
        - values : 'symb', or 'num'.
            * If  'sym', return the symbolic coefficients.
            * If  'num', return the numerical coefficients for the default values of the components set in the Circuit object.
        - scheme (str): Discretization scheme ('frwrd', 'bckwrd', 'blnr').
        - Fs (float or sympy.Symbol): Samplerate.
        - norm (bool): Returns normalized coefficients w.r.t a[0].
        - simple (bool): Analog coefficients are expressed implicitly, returns generic discretized coefficients.

        Returns:
        - Bd (np.array): Discretized numerator coefficients.
        - Ad (np.array): Discretized denominator coefficients.
        """
        if values == 'sym':
            Bd, Ad = self.sym_coeffz(scheme=scheme, norm=norm, simple=simple)
        elif values == 'num':
            Bd, Ad = self.num_coeffz(scheme=scheme, Fs=Fs)
        else:
            raise ValueError('values has to be set to sym or num')
            
        return Bd, Ad

    def __str__(self):
        return str(self.sympyExpr)
    
    def __repr__(self):
        return repr(self.sympyExpr)
    
    def __getattr__(self, name):
        '''
        I want to be able to still use the sympy library to manipulate the symbolic transfer function, 
        that's why I override the __getattr__ method to get the attribute of the symbolic transfer function directly
        '''
        return getattr(self.sympyExpr, name)
    
    def __call__(self, *args, **kwargs):
        return self.sympyExpr(*args, **kwargs)



def plotTransfertFunction(f, h, legend=None, title=None, semilogx=True, dB=True, phase=True):
    ''' 
    Plot the frequency response of the filter.

    Parameters:
    - f (array): frequency array (Hz) -> 1D array
    - h (array): frequency response array (complex) -> Max 3D array!
    - legend (dict or str): Legend of the different lines.
                            Use the same argument as component_values in the 
                            num_coeffs() method if 
                            you want a fast legend creation!
    - title (str): title of the plot
    - semilogx (bool): If True, use a logarithmic scale for the x-axis
    - semilogy (bool): If True, use a logarithmic scale for the y-axis
    - dB (bool): If True, plot the magnitude in dB
    - phase (bool): If True, plot the phase in radians
    '''
    
    if h.ndim > 3:
        raise ValueError("The input array h must have at most 3 dimensions.")

    #default colors and linestyles 
    colors = plt.cm.tab10(np.linspace(0, 1, 10)) #RdYlBu, plasma, viridis, inferno, magma, cividis
    linestyles = ['-', '--', '-.', ':', '-*', '-o', '-+', '-x', '-|']
    colors_nbr = len(colors)
    linestyles_nbr = len(linestyles)

    h = np.atleast_2d(h)

    if dB: 
        h_mag =  20 * np.log10(np.abs(h))
        scale = 'dB'
    else:
        h_mag = np.abs(h)
        scale = 'lin.'

    fig, ax = plt.subplots(2 if phase else 1, 1, sharex=True, layout='tight')
    ax = np.atleast_1d(ax)

    if isinstance(title, str):
        ax[0].set_title(title)

    if semilogx:
        ax[0].set_xscale('log')

    plot_curves(ax[0], f, h_mag, linestyles, colors)

    ax[0].set_ylabel(f'Magnitude [{scale}]')
    ax[0].grid(which='both')
    ax[0].format_coord = lambda x, y: f'x={x:.3f}, y={y:.3f}'

    if phase:
        
        plot_curves(ax[1], f, np.angle(h), linestyles, colors)
        
        ax[1].set_ylabel('Phase [radians]')
        ax[1].set_xlabel('Frequency [Hz]')
        ax[1].set_yticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
        ax[1].set_yticklabels(['$-\pi$', '$-\pi/2$', '0', '$\pi/2$', '$\pi$'])
        ax[1].set_xlim([f[0], f[-1]])
        ax[1].format_coord = lambda x, y: f'x={x:.3f}, y={y:.2f}'
        ax[1].grid(which='both')
    
    ax[-1].set_xticks([1, 20, 50, 100, 200, 500, 1e3, 2e3, 5e3, 10e3, 20e3, 100e3, 1e6])
    ax[-1].set_xticklabels(['1', '20', '50', '100', '200', '500', '1k', '2k', '5k', '10k', '20k', '100k', '1M'])
    ax[-1].set_xlim([f[0], f[-1]])

    plot_legend(ax, h, legend, colors, linestyles)
    plt.show()


def plot_curves(ax, x, y, linestyles, colors):
    if y.ndim == 2:
        for j, y_array in enumerate(y):
            ax.plot(x, y_array, color=colors[j % len(colors)])
    else:
        for i, y_slice in enumerate(y):
            for j, y_array in enumerate(y_slice):
                ax.plot(x, y_array, linestyles[i % len(linestyles)], color=colors[j % len(colors)])
    
    

def plot_legend(ax, h, legend, colors, linestyles):

    colors_nbr, linestyles_nbr = len(colors), len(linestyles)

    if legend is not None: 

        if isinstance(legend, str):
            ax[0].plot([], [], '-', label=legend)
            
        else:
            if h.ndim == 3:
                first_key = list(legend.keys())[0]
                first_key_values = legend[first_key]

                for i in range(len(first_key_values)):
                    ax[0].plot([], [], linestyles[i % linestyles_nbr], color='grey', label=f'{first_key}: {first_key_values[i]}')
    
            last_key = list(legend.keys())[-1]  #if h.ndim==2, last_key is obviously the only key in the legend dictionary
            last_key_values = legend[last_key]

            for i in range(len(last_key_values)):
                ax[0].plot([], [], '-', color=colors[i % colors_nbr], label=f'{last_key}: {last_key_values[i]}')
        ax[0].legend(loc='best')
