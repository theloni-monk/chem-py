import sys
import numpy as np
import sympy as sm

_DEBUG = False



def dprint(*x):
    if _DEBUG: print(*x)


def balanceEq(EQ, getVect = False):
    eq = EQ
    bEQ = ""
    dprint(eq)

    reactants, products = eq.split('->')[0], eq.split('->')[1]
    reactants = str(reactants).replace(' ', '').split('+')
    products = str(products).replace(' ', '').split('+')
    dprint(reactants, products)

    # double check user input if there are the same elements on each side
    assertstring = "Equation cannot be balanced: " + eq + "\n Reactants: " + str(reactants) + "\t Products: " + str(
        products) + "\n Reactants: " + str(getElements(reactants)) + "\t Products: " + str(getElements(products))
    assert getElements(reactants) == getElements(products), assertstring

    elements = getElements(reactants)
    dprint(elements)

    # compute matrix
    compound_matrix = []
    for compound in reactants:
        compound_matrix.append(getElementVector(compound, elements))

    for compound in products:
        compound_matrix.append([-i for i in getElementVector(compound, elements)])
    
    # compute reduced row echalon form
    compound_matrix = sm.Matrix(compound_matrix).T
    #sm.pprint(compound_matrix)
    rrMat = compound_matrix.rref()[0]
    #sm.pprint(rrMat)
    dprint(rrMat.shape)

    # substitute lcm of values
    # flatten solution matrix, add last term for the reduced form
    solution_vect = [-i[0] for i in rrMat.col(rrMat.shape[1] - 1).tolist()]
    dprint(solution_vect)

    # shitty inefficient hack because im lazy
    mults = []
    for i in solution_vect:
        for j in range(1, 50):
            if int(j*i) == j*i:  # lazy check for integer lcm
                mults.append(j)
                break

    lcm = np.lcm.reduce(mults)
    solution_vect = [int(lcm*i) for i in solution_vect]
    solution_vect.append(lcm)

    if getVect: return solution_vect
    
    # reconstruct balanced equation
    # TODO: print prettier
    # BUG: in case that initial coefficient is given, it will write 23Fe instead of 6Fe or 2*3Fe given initial compound 3Fe
    for c in range(len(reactants)):
        bEQ += (str(solution_vect[c]) if solution_vect[c] != 1 else '') + (reactants[c]) + " + "
    bEQ = bEQ[:-3] + ' -> '
    for c in range(len(reactants), len(reactants)+len(products)):
        bEQ += (str(solution_vect[c]) if solution_vect[c] != 1 else '' )+ (products[c-len(reactants)]) + " + "
    dprint(solution_vect)
    return bEQ[:-3]  # -3 removes the extra ' +  '


# HELPERS:
def getElements(compounds):  # get elements and sort by first char value of element
    elements = []
    for compound in compounds:
        for c in range(len(compound)):
            # detect capital letter: start of element
            if ord(compound[c]) > 64 and ord(compound[c]) < 91:
                try:
                    # if capital has lowercase after it concat into element
                    if ord(compound[c+1]) > 96:
                        elements.append(compound[c]+compound[c+1])
                    else:
                        # otherwise just append capital letter to element list
                        elements.append(compound[c])
                except IndexError:
                    # catch single-letter element at end
                    elements.append(compound[c])
                    break  # if you reach the end of the string
    elements.sort()  # sort elements so sets are always same and can be compared

    # purge copies of same element across compounds
    i = 0
    while i < len(elements):
        try:
            if elements[i] == elements[i+1]:
                elements.remove(elements[i])
            else:
                i += 1
        except IndexError:
            break

    return elements


def getLeadingCof(string):
    cof = ""

    if not len(string): return 1

    try:
        for i in range(len(string)):  # append number until you hit an invalid int literal
            # throws value error if string[i] not int literal
            x = int(string[i])
            cof += string[i]  # concat string
    except ValueError:
        try:
            cof = int(cof)
        except ValueError:
            cof = 1  # if there weren't any cofs then set it to 1
    return int(cof)


def getElementVector(c, elements):
    dprint('getElementVector: ')
    dprint(c)
    compound = c
    compoundEls = getElements([compound])
    compound_vect = [0 for i in elements] # zeros vector of elements

    compoundCof = getLeadingCof(compound)  # get molecular cof
    if compoundCof != 1: compound = compound[len(str(compoundCof)):]  # slice off leading cof if there is one
    # count elements to form vector
    for i in range(len(elements)):
        dprint('per element')
        dprint(elements[i], elements, compound)
        if elements[i] in compoundEls:  # iterate through elements in compound
            # the cof after each element in the compound
            compound = compound.split(elements[i])[1]
            dprint(compound)
            compound_vect[i] = getLeadingCof(compound)
            if compound_vect[i] != 1: compound = compound[len(str(compound_vect[i])):]
            dprint(compound_vect[i])

    dprint(compound_vect)

    return [ compoundCof * i for i in compound_vect]

# MAIN
if __name__ == "__main__":
    if len(sys.argv)>1:
        print('\tSolved EQ: | ' + balanceEq(sys.argv[1]) + ' |')
    else:
        print(balanceEq(input("Input equation: ")))
