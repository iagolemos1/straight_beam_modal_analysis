import numpy as np
import matplotlib.pyplot as plt

class straight_beam(object):
    """
    This class refers to the characterisation of a straight beam, 
    with no section variation, for an analytical analysis and study pourposes.
    
	Parameters
	----------
    mat_YoungMod : float
        Young's modulus of the material (default: Structural steel - 2.1e11 Pa)
    mat_density : float
        Density of the material (default: Structural steel - 7800 kg/m^3)
    mat_poisson : float
        Poisson ratio of the material (default: Structural steel - 0.3)
	beam_length : float
		Length of the beam (default = 1 m)
    mom_in : float
        Moment of inertia of the section (default : None)
    sec_area : float
        Cross-section area of the beam (default : None)
        
	Public attributes
	-----------------
	checked : bool
		Status of validation of arguments. After validation, if parameters obbey their restrictions, it will be 'True'. Else, it will remain 'False'.
	E : float
		Variable associated to the Young's modulus of the material
	rho : float
		Variable associated to the density of the material
	nu : float
		Variable associated to the poisson ratio of the material
    G : float
        Shear modulus of the material
	L : float
		Variable associated to the length of the beam
    I : float
        Moment of inertia of the section, used for bending calculation
    A : float
        Cross-section area of the beam, used for beding calculation
    """
    
    def __init__(self, mat_YoungMod = 2.1e11, mat_density = 7800, mat_poisson = 0.3, beam_length = 1, mom_in = None):
        self.checked = False
        self.E = mat_YoungMod
        self.rho = mat_density 
        self.nu = mat_poisson
        self.L = beam_length
        self.G = self.E/(2*(1 + self.nu))
        self.I = mom_in
        self.A = sec_area
        
    def verifyArgs(self):
        if self.E and self.rho and self.nu and self.L > 0:
            self.checked = True
            
    def knL_long_tors(self, n, limite_0, limite_1):
        """
        Method to compute the wave number
        ----------
        n : list or nummpy.ndarray
            Modes to compute
        limite_0 : str
            Structure constraint in x = 0. Can be 'fix' or 'free'. 
        limite_1 : str
            Structure constraint in x = L. Can be 'fix' or 'free'. 
        Returns
        -------
        knL : numpy.ndarray
            Wave number multiplied by the length of the beam.
        """
        if type(n).__module__ != np.__name__:
            n = np.array(n)
        if limite_0 == 'fix' and limite_1 == 'fix':
            knL = n*np.pi
        elif limite_0 == 'fix' and limite_1 == 'free' or limite_0 == 'free' and limite_1 == 'fix':
            knL = (n - 0.5)*np.pi
        elif limite_0 == 'free' and limite_1 == 'free':
            knL = np.concatenate((np.array([0]), n*np.pi))
            print(knL)
        return(knL)
    
    def knL_flex(self, n, limite_0, limite_1):
        """
        Method to compute the wave number
        ----------
        n : list or nummpy.ndarray
            Modes to compute
        limite_0 : str
            Structure constraint in x = 0. Can be 'fix', 'free' or 'supported'. 
        limite_1 : str
            Structure constraint in x = L. Can be 'fix', 'free' or 'supported'. 
        Returns
        -------
        knL : numpy.ndarray
            Wave number multiplied by the length of the beam.
        """
        if type(n).__module__ != np.__name__:
            n = np.array(n)
        if limite_0 == 'fix' and limite_1 == 'fix':
            knL = (2*n + 1)*np.pi/2
        elif limite_0 == 'fix' and limite_1 == 'free' or limite_0 == 'free' and limite_1 == 'fix':
            knL = (2*n - 1)*np.pi/2
        elif limite_0 == 'fix' and limite_1 == 'supported' or limite_0 == 'supported' and limite_1 == 'fix':
            knL = (4*n + 1)*np.pi/4
        elif limite_0 == 'free' and limite_1 == 'free':
            knL = np.concatenate((np.array([0]), (2*n + 1)*np.pi/2))
        elif limite_0 == 'free' and limite_1 == 'supported' or limite_0 == 'supported' and limite_1 == 'free':
            knL = np.concatenate((np.array([0]), (4*n + 1)*np.pi/4))
        elif limite_0 == 'supported' and limite_1 == 'supported':
            knL = n*np.pi
        return(knL)
    
    def modal_traction(self, n = [1,2,3], limite_0 = 'fix', limite_1 = 'free', plot_defs = True):
        """
        Method to compute the modal analysis of the beam in the case of traction
        ----------
        n : list or nummpy.ndarray
            Modes to compute. The default is [1, 2, 3] (computes the first, second and third mode).
        limite_0 : str
            Structure constraint in x = 0. Can be 'fix' or 'free'. The default is 'fix'. 
        limite_1 : str
            Structure constraint in x = L. Can be 'fix' or 'free'. The default is 'free'.
        plot_defs: booleon, optional
            Wheter to plot or not the mode shapes.
        Returns
        -------
        eigen_freqs : numpy.ndarray
            Eigenfrequencies in Hz of the given structure.
        """
        c = (self.E/self.rho)**0.5
        knL = self.knL_long_tors(n, limite_0, limite_1)
        eigen_freqs = (knL/(self.L))*c/(2*np.pi)
        if plot_defs == True:
            x = np.linspace(0, self.L, 500)
            deformations = np.zeros((np.shape(eigen_freqs)[0], np.shape(x)[0]))
            for index, val in enumerate(eigen_freqs):
                if limite_0 == 'fix' and limite_1 == 'fix':
                    deformations[index, :] = np.sin((knL/(self.L))[index]*x)
                elif limite_0 == 'fix' and limite_1 == 'free' or limite_0 == 'free' and limite_1 == 'fix':
                    deformations[index, :] = np.sin((knL/(self.L))[index]*x)
                elif limite_0 == 'free' and limite_1 == 'free':
                    deformations[index, :] = np.cos((knL/(self.L))[index]*x)
            plt.figure()
            for i in range(len(n)):
                plt.plot(x, deformations[i,:], label = ('$\u03C6_{}$ : {:.2f} Hz'.format(i + 1, eigen_freqs[i])))
                plt.legend()
                plt.ylabel('$\u03C6_{n}(x)$')
                plt.xlabel('x')
                plt.xlim(0,self.L)
            plt.title('Traction mode shapes')
        return(eigen_freqs)
    
    def modal_torsion(self, n = [1,2,3], limite_0 = 'fix', limite_1 = 'free', plot_defs = True):
        """
        Method to compute the modal analysis of the beam in the case of traction
        ----------
        n : list or nummpy.ndarray
            Modes to compute. The default is [1, 2, 3] (computes the first, second and third mode).
        limite_0 : str
            Structure constraint in x = 0. Can be 'fix' or 'free'. The default is 'fix'. 
        limite_1 : str
            Structure constraint in x = L. Can be 'fix' or 'free'. The default is 'free'.
        plot_defs: booleon, optional
            Wheter to plot or not the mode shapes.
        Returns
        -------
        eigen_freqs : numpy.ndarray
            Eigenfrequencies in Hz of the given structure.
        """
        c = (self.G/self.rho)**0.5
        knL = self.knL_long_tors(n, limite_0, limite_1)
        eigen_freqs = (knL/(self.L))*c/(2*np.pi)
        if plot_defs == True:
            x = np.linspace(0, self.L, 500)
            deformations = np.zeros((np.shape(eigen_freqs)[0], np.shape(x)[0]))
            for index, val in enumerate(eigen_freqs):
                print(deformations[index,:])
                if limite_0 == 'fix' and limite_1 == 'fix':
                    deformations[index, :] = np.sin((knL/(self.L))[index]*x)
                elif limite_0 == 'fix' and limite_1 == 'free' or limite_0 == 'free' and limite_1 == 'fix':
                    deformations[index, :] = np.sin((knL/(self.L))[index]*x)
                elif limite_0 == 'free' and limite_1 == 'free':
                    deformations[index, :] = np.cos((knL/(self.L))[index]*x)
            plt.figure()
            for i in range(len(n)):
                plt.plot(x, deformations[i,:], label = ('$\u03C6_{}$ : {:.2f} Hz'.format(i + 1, eigen_freqs[i])))
                plt.legend()
                plt.ylabel('$\u03C6_{n}(x)$')
                plt.xlabel('x')
                plt.xlim(0,self.L)
            plt.title('Torsion mode shapes')
        return(eigen_freqs)
    
    def modal_blend(self, n = [1,2,3], limite_0 = 'fix', limite_1 = 'free', plot_defs = True):
        """
        Method to compute the modal analysis of the beam in the case of bending
        ----------
        n : list or nummpy.ndarray
            Modes to compute. The default is [1, 2, 3] (computes the first, second and third mode).
        limite_0 : str
            Structure constraint in x = 0. Can be 'fix' or 'free'. The default is 'fix'. 
        limite_1 : str
            Structure constraint in x = L. Can be 'fix' or 'free'. The default is 'free'.
        plot_defs: booleon, optional
            Wheter to plot or not the mode shapes.
        Returns
        -------
        eigen_freqs : numpy.ndarray
            Eigenfrequencies in Hz of the given structure.
        """
        c = (self.E*self.I/(self.rho*self.A))**0.5
        kn = knL/self.L
        knL = self.knL_flex(n, limite_0, limite_1)
        eigen_freqs = (c**0.5)*(kn**2)
        if plot_defs == True:
            x = np.linspace(0, self.L, 500)
            deformations = np.zeros((np.shape(eigen_freqs)[0], np.shape(x)[0]))
            for index, val in enumerate(eigen_freqs):
                print(deformations[index,:])
                if limite_0 == 'fix' and limite_1 == 'fix':
                    deformations[index, :] = np.sin(kn[index]*x)
                elif limite_0 == 'fix' and limite_1 == 'free' or limite_0 == 'free' and limite_1 == 'fix':
                
                elif limite_0 == 'fix' and limite_1 == 'supported' or limite_0 == 'supported' and limite_1 == 'fix':
                
                elif limite_0 == 'free' and limite_1 == 'free':
                
                elif limite_0 == 'free' and limite_1 == 'supported' or limite_0 == 'supported' and limite_1 == 'free':
                
                elif limite_0 == 'supported' and limite_1 == 'supported':
                    
            plt.figure()
            for i in range(len(n)):
                plt.plot(x, deformations[i,:], label = ('$\u03C6_{}$ : {:.2f} Hz'.format(i + 1, eigen_freqs[i])))
                plt.legend()
                plt.ylabel('$\u03C6_{n}(x)$')
                plt.xlabel('x')
                plt.xlim(0,self.L)
            plt.title('Torsion mode shapes')
        return(eigen_freqs)   

        
        
    
    
    