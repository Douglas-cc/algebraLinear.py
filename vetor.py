import numpy as np
from math import sqrt, acos, pi
from decimal import Decimal, getcontext

getcontext().prec = 5

class Vetor(object):

    CANNOT_NORMALIZE_ZERO_VECTOR_MSG = 'Não é possível normalizar o vetor nulo'

    def __init__(self, coordenadas):
        try:
            if not coordenadas:
                raise ValueError

            self.coordenadas = list(Decimal(x) for x in coordenadas)             
            self.dimensao = len(self.coordenadas)

        except ValueError:
            raise ValueError('As coordenadas devem ser vazias')

        except TypeError:   
            raise TypeError('As coordenadas devem ser vazias')                    


    def soma_vt(self, vetor):
        resultado = [x + y for x,y in zip(self.coordenadas, vetor.coordenadas)]
        return Vetor(resultado)        


    def subt_vt(self, vetor):
        resultado = [x - y for x,y in zip(self.coordenadas, vetor.coordenadas)]
        return Vetor(resultado)


    def pdt_escalar_vt(self, vetor): # conhecido também como produto interno do vetor
        return sum([x*y for x,y in zip(self.coordenadas, vetor.coordenadas)])


    def mult_escalar(self, escalar): 
        resultado = [Decimal(escalar) * x for x in self.coordenadas]
        return Vetor(resultado)
  
    def vt_unitario(self):
        if len(self.coordenadas) == 1:
            return f'Vector: {self.coordenadas} é um vetor unitario'
        elif sum(self.coordenadas) == 1:
            return f'Vector: {self.coordenadas} é um vetor unitario'
        else:    
            return f'Vector: {self.coordenadas} não é um vetor unitario'

    def ponto_medio_vt(self, vetor):
        resultado = [(x + y)/2 for x,y in zip(self.coordenadas, vetor.coordenadas)]
        return Vetor(resultado)

    # Basicamente aplicamos o teorema de pitagoras nos valores do array para obter magnitude
    def magnitude_vt(self):
        coord_ao_quadrado = [x**2 for x in self.coordenadas]
        return Decimal(sqrt(sum(coord_ao_quadrado)))

    #  Primeiro garantimos se o produto escalar esteja enre 10 e -10, para em seguida
    #  descubro se a magnitude do vetor é zero, se sim é o vetor é nulo
    def vt_zero(self, tolerancia=1e-10):
        return self.magnitude_vt() < tolerancia

    '''
    Normalizar é o processo de encontrar um vetor unitarios na mesma direção que um vetor doador
    Para isso encontramos sua magnitude e depois multplicamos por um escalar ou seja multplicamos
    o vetor por 1 sobre a magnitude do vetor, isso dimensiona o vetor para que ele seja iqual a 1
    '''
    def normaliza_vt(self):
        try:
            return self.mult_escalar(Decimal('1.0')/self.magnitude_vt())
        except ZeroDivisionError:
            raise Exception(CANNOT_NORMALIZE_ZERO_VECTOR_MSG)        

    def angulo_vt(self, vetor, em_graus=True):
        try:
            vetor1 = self.normaliza_vt()
            vetor2 = vetor.normaliza_vt()
            radianos = acos(vetor1.pdt_escalar_vt(vetor2)) # angulo em radianos
            
            # usamos essa variavel para auxiliar na conversão de radianos para graus
            if em_graus:
                return radianos * 180/pi
            else:
                radianos

        except Exception as e:
            if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
                raise Exception('Não é possível calcular um ângulo com o vetor nulo')
            else:
                raise e

    #  Primeiro garantimos se o produto escalar esteja enre 10 e -10, para em seguida
    #  e por fim descubro se o produto interno deles é zero, se sim são ortogonais
    def vt_ortagonal(self, vetor, tolerancia=1e-10):
        return abs(self.pdt_escalar_vt(vetor)) < tolerancia

    # Dado algumas condições para determinar se meu vetor é paralelo,se os vetores
    # são zero, se o angulo é zero ou se o angulo tem o valor de pi radianos
    def vt_paralelo(self, vetor):
        return(self.vt_zero() or
                vetor.vt_zero() or
                self.angulo_vt(vetor) == 0 or
                self.angulo_vt(vetor) == pi)  

    '''
    Calculamos a projeção de um vetor sobre o vetor base normalizando o vetor base 
    para formar um vetor unitario, em seguida calcula o produto escalar do vetor 
    unitario e multiplico o vetor unitario pelo seu produto escalar produto escalar  
    '''
    def comp_paralelo(self, vetor_base):    
        try:
            vetor_unitario = vetor_base.normaliza_vt()
            produto_escalar = self.pdt_escalar_vt(vetor_unitario)
            return vetor_unitario.multi_escalar(produto_escalar)
        # Se vetor da base for zero,retornara esta 
        # exceção ou seja normalização do vetor falhou
        except Exception as e:
            if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
                raise Exception(self.NO_UNIQUE_PARALLEL_COMPONENT_MSG)
            else:
                raise e                

    # Para calcular o componente ortogonal é preciso calcula a projeção,
    # que é component paralelo do vetor, base menos o vetor base
    def comp_ortogonal(self, vetor_base):
        try:
            projecao = self.comp_paralelo(vetor_base)
            return self.subt_vt(projecao)
        except Exception as e:
            if str(e) == self.NO_UNIQUE_PARALLEL_COMPONENT_MSG:
                raise Exception(self.NO_UNIQUE_ORTHOGONAL_COMPONENT_MSG)
            else:
                raise e    
    
    # prduto vetorial, para obetermos o novo vetor praticamente calculamos uma determinante
    def pdt_vetorial(self, vetor):
        try:
            x1, y1, z1 = self.coordenadas
            x2, y2, z2 = vetor.coordenadas
            novo_vetor = [ y1*z2 - y2*z1,
                        -(x1*z2 - x2*z1),
                        x1*y2 - x2*y1 ]

            return Vetor(novo_vetor)             

        except ValueError as e:
            msg = str(e)
            if msg == 'need more than 2 values to unpack':
                self_em_R3 = Vetor(self.coordenadas + ('0',))
                vetor_em_R3 = Vetor(vetor.coordenadas + ('0',))
                return self_em_R3.pdt_vetorial(vetor_em_R3)
            elif (msg == 'too many values to unpack' or
                  msg == 'need more than 2 values to unpack'):
                raise Exception(self.ONLY_DEFINED_IN_TWO_THREE_DIMS_MSG)      
            else:
                raise e

    '''
    O produto Kronecker, denotado por ⊗, é uma operação em duas matrizes de tamanho arbitrário
    resultando em uma matriz de bloco. É uma generalização do produto externo  de vetores para 
    matrizes, fornece a matriz do produto tensorial com relação a uma escolha padrão de base
    '''
    def pdt_kronecker(self, vetor):
        return np.kron(self.coordenadas, vetor.coordenadas)

    # area do paralelograma é magnitude do produto vetorial do vetor
    def area_paralelograma(self, vetor):
        return pdt_vetorial(vetor).magnitude_vt()

    # area do triagulo formado pelos vetore é area do paralelograma  dividido por dois
    def area_triagulo(self, vetor):
        return self.area_paralelograma(vetor)/Decimal('2.0')

    # Quando escrevemos uma nova classe, usamos o __init__,que facilita 
    # a instanciação de objetos e __str__, que é útil para a depuração.            
    def __str__(self):
        return f'Vetor: {self.coordenadas}'

    # Aqui estmos comparando os vetores instanciados, para ver se são iguais
    # ou seja se possuem o mesmo valor de deslocamento em cada eixo
    def __eq__(self, vetor):
        return self.coordenadas == vetor.coordenadas