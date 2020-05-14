from math import sqrt, acos, pi 
from decimal import Decimal, getcontext

getcontext().prec = 30

class Vector(object):

    CANNOT_NORMALIZE_ZERO_VECTOR_MSG = 'Cannot normalize the zero vector'
    
    def __init__(self, coordenadas):
        try:
            if not coordenadas:
                raise ValueError
            self.coordenadas = tuple([Decimal(x) for x in coordenadas])
            self.dimensao = len(self.coordenadas)

        except ValueError:
            raise ValueError('As coordenadas devem ser vazias')

        except TypeError:   
            raise TypeError('As coordenadas devem ser vazias')
        
    def soma(self, v):
        nova_coordenada = [x+y for x,y in zip(self.coordenadas, v.coordenadas)]
        return Vector(nova_coordenada)        
       
    def subtracao(self, v):
        nova_coordenada = [x-y for x,y in zip(self.coordenadas, v.coordenadas)]
        return Vector(nova_coordenada)    

    def multi_escalar(self, c):    
        nova_coordenada = [Decimal(c)*x for x in self.coordenadas]
        return Vector(nova_coordenada)
    
    def magnitude(self):
        coord_quadrado = [x**2 for x in self.coordenadas]
        return Decimal(sqrt(sum(coord_quadrado)))
    
    def normalizacao(self):
        try:
            magnitude = self.magnitude()
            return self.multi_escalar(Decimal('1.0')/self.magnitude())
        except ZeroDivisionError:
            raise Exception(CANNOT_NORMALIZE_ZERO_VECTOR_MSG)
    
    def produto_interno(self, v):
        return sum([x*y for x,y in zip(self.coordenadas, v.coordenadas)])
    
    def angulo(self, v, graus=False):
        try:
            v1 = self.normalizacao()
            v2 = v.normalizacao()
            radianos = acos(v1.produto_interno(v2))
           
            if graus:
                converter = 180 / pi
                return radianos * converter
            else:
                return radianos

        except Exception as e:
            if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
                raise Exception('Não é possível calcular um ângulo com o vetor zero')
            else:
                raise e

    def vt_ortagonal(self,v, tolerancia=1e-10):
        return abs(self.produto_interno(v) < tolerancia)

    def vt_paralelo(self,v):
        return (self.vt_zero() or 
                v.vt_zero() or 
                self.angulo(v) == 0 or 
                self.angulo(v) == pi )  
    
    def vt_zero(self, tolerancia=1e-10):
        return self.magnitude() < tolerancia
    
    def componente_ortagonal(self, base):
        try:
            projecao = self.componente_paralelo(base)
            return self.subtracao(projecao)
        
        except  Exception as e:
            if str(e) == self.NO_UNIQUE_PARALLEL_COMPONENT_MSG:
                raise Exception(self.NO_UNIQUE_ORTHOGONAL_COMPONENT_MSG)
            else:
                raise e
    def componente_paralelo(self, base):    
        try:
            u = base.normalizacao()
            peso = self.produto_interno(u)
            return u.multi_escalar(peso)

        except Exception as e:
            if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
                raise Exception(self.NO_UNIQUE_PARALLEL_COMPONENT_MSG)
            else:
                raise e

    def __str__(self):
        return 'Vector: {}'.format(self.coordenadas)

    def __eq__(self, v):
        return self.coordenadas == v.coordenadas

print('\t#3')
v = Vector([3.009, -6.172, 3.692, -2.51])
w = Vector([6.404, -9.144, 2.759, 8.718])

v.componente_paralelo(w)
