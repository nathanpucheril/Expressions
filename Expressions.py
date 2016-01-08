import utils
import math
from copy import deepcopy

# @author: Nathan Pucheril

class ExpressionBuilder(object):

    @staticmethod
    def parse(string):
        pass

    @staticmethod
    def make_polynomial(terms):
        return Polynomial(terms)

    def __div__(num, den):
        # handle cases of inserting constants and expressions
        return Fraction(num, den)

    @staticmethod
    def x(p):
        return PowerTerm(1, X(), p)

    @staticmethod
    def log(x):
        # if not Expression.isExpression(x):
        #     Expressify

        return LogTerm(x)

    @staticmethod
    def ln(x):
        return LogTerm(x, 1, math.e)

    @staticmethod
    def log_base(base, x):
        return LogTerm(x, 1, base)

    def __add__(x, y):
        return x + y

    def __mul__(x, y):
        return x * y

    @staticmethod
    def exp(x):
        return ExponentialTerm(1, x)

    @staticmethod
    def cons(constant):
        return constant

    @staticmethod
    def coefficient(coefficient, term):
        c = Expression.make_constant(constant)
        term.set_coefficient(c)
        return term

#############################

class Expression(object):

    def __init__(self, expressions):
        assert isinstance(expressions, list), "Terms must be a list"
        assert all(map(Expression.isExpression, expressions)
                   ), "All Terms in an Expression must be of the Class Expression"
        self.termsList = expressions
        self._termFlatten()
        self.c = 1

    def _termFlatten(self):
        terms = []
        repeat = 0
        for term in self.termsList:
            if term.__class__ == Expression:
                terms.extend(term.terms)
                repeat = 1
            else:
                terms.append(term)
        self.termsList = terms
        if repeat:
            self._termFlatten()

    @property
    def terms(self):
        return self.termsList

    @property
    def coefficient(self):
        return self.c

    def copy(self):
        return deepcopy(self)

    def set_coefficient(self, coefficient):
        self.c = coefficient
        return self

    def _clean(self):
        cleaned = []
        for term in self.terms:
            c = term.coefficient
            if c == 0:
                continue
            cleaned.append(term)
        self.terms = cleaned
        return self

    def simplify(self):
        return self._simplifier()
        # return self

    def _simplifier(self):
        self.termsList = [term if utils.isnumeric(term) else term.simplify() for term in self.terms]
        self.termsList = list(filter(lambda x: x, self.terms)) ## NEED TO FIX
        return self

    def __call__(self, x):
        assert utils.isnumeric(x)
        return self.evaluate(x)

    def __add__(e1, e2):
        assert map(Expression.isExpression, (e1, e2)
                   ), "Arguements must be an Expression"
        return Expression([e1 , e2])

    def __radd__(e1,e2):
        return e1.__add__(e2)

    def __sub__(e1, e2):
        # ERROR HANDLING DONE IN ADD
        return e1 + -e2

    def __rsub__(e1,e2):
        return e1.__sub__(e2)

    def __mul__(e1, e2):
        assert map(Expression.isExpression, (e1, e2)
                   ), "Arguements must be of Class Term"
        multiplied = 0
        for term1 in [e1] if utils.isnumeric(e1) else e1.terms:
            for term2 in [e2] if utils.isnumeric(e2) else e2.terms:
                multiplied += (term1 * term2)
        return multiplied

    def __rmul__(e1, e2):
        return e1.__mul__(e2)

    def __div__(e1, e2):
        """ p1  divided by p2 -> p1/p2 """
        # Error Handling done in Mul
        if utils.isnumeric(e2):
            return e1 * (1.0/e2)
        return e1 * e2.reciprocal()

    def __rdiv__(e1,e2):
        return e1.__div__(e2)

    def __neg__(self):
        # ERROR HANDLING DONE IN MUL
        return self.copy() * -1

    def __pow__(e1, e2):
        return PowerTerm(1, e1, e2)

    def __rpow__(e1, e2):
        return e1.__pow__(e2)

    def reciprocal(self):
        return Fraction._fractifyExpression(self).invert()

    def invert(self):
        # Not In place because if it isnt of type fraction, must return a
        # fraction
        return self.reciprocal()

    def derivative(self):
        return Expression([term.derivative for term in self.terms])

    @staticmethod
    def isZero(self):
        return False

    @staticmethod
    def isExpression(e):
        return isinstance(e, Expression) or utils.isnumeric(e)

    def evaluate(self, x):
        return sum([term(x) for term in self.termsList])

    def toPolynomial(self):
        return Polynomial(self.terms)

    def __str__(self):
        output = str(self.c) if self.c != 1 else ""
        for term in self.terms:
            output += str(term) + " + "
        return output[:len(output) - 3]


class Fraction(Expression):

    def __init__(self, num, den, ):
        assert map(Expression.isExpression, (num, den))
        assert not Expression.isZero(den)
        # CHECK FOR SIMPLIFIYING FRACTIONS
        self.num = num
        self.den = den
        self.reduceCoefficients()

    @property
    def terms(self):
        return [self]

    @staticmethod
    def _fractifyExpression(e):
        if Fraction.isFraction(e):
            return e
        return Fraction(e, 1)

    def simplify(self):
        return self

    def __add__(f1, f2):
        assert isinstance(f1, Expression) and isinstance(f2, Expression)
        f1 = Fraction._fractifyExpression(f1)
        f2 = Fraction._fractifyExpression(f2)
        return Fraction(f1.num * f2.den + f2.num * f1.den, f1.den * f2.den)

    def __mul__(f1, f2):
        assert isinstance(f1, Expression) and isinstance(f2, Expression)
        f1 = Fraction._fractifyExpression(f1)
        f2 = Fraction._fractifyExpression(f2)
        return Fraction(f1.num * f2.num, f1.den * f2.den)

    def __div__(f1, f2):
        assert isinstance(f1, Fraction) and isinstance(f2, Fraction)
        return f1 * f2.copy().invert()
    #
    # def __neg__(f1):
    #     return Fraction(-self.num, self.den)

    def invert(self):
        self.num, self.den = self.den, self.num
        return self

    def reciprocal(self):
        return Fraction(self.den, self.num)

    def reduceCoefficients(self):
        coefficients = list(
            term.coefficient for term in self.num.terms + self.den.terms)

        def gcd_help(x, y):
            x, y = map(abs, (x, y))
            if y > x:
                return gcd_help(y, x)
            if y == 0:
                return x
            return gcd_help(y, x % y)
        gcd = reduce(gcd_help, coefficients)
        for x in self.num.terms:
            x.set_coefficient(x.coefficient / gcd)
        for x in self.den.terms:
            x.set_coefficient(x.coefficient / gcd)
        return self

    def evaluate(self, x):
        return self.num(x) / self.den(x)

    @staticmethod
    def isFraction(f):
        return isinstance(f, Fraction)

    def __str__(self):
        if str(self.den) == "1.0" or str(self.den) == "1":
            return str(self.num)
        return "(" + str(self.num) + ")/(" + str(self.den) + ")"


class Polynomial(Expression):

    def __init__(self, terms):
        assert isinstance(terms, list), "Terms must be a list of tuples"
        for term in terms:
            assert isinstance(
                term, PowerTerm), "Every Term must be a PowerTerm"
            assert utils.isnumeric(term.exp), "Every Term must be a PowerTerm"
        super(Polynomial, self).__init__(terms)
        self.termsList = terms

    def simplify(self):
        return self._combine_terms()._simplifier()._sort()

    def _sort(self):
        self.termsList = sorted(self.terms, key=lambda x: x.exp, reverse=True)
        return self

    def _combine_terms(self):

        return self

    def _simplifier(self):

        self.termsList = list(filter(lambda x: x, self.terms))
        return self

    # Done
    def __str__(self):
        output = ""
        for term in self.terms:
            c, exp = term.c, term.exp
            if c == 0:
                continue
            elif exp == 0:
                output += str(c) + " + "
            elif c == 1:
                output +="^" + str(exp) + " + "
            else:
                output += str(c) +  "^" + str(exp) + " + "
        output = output[: len(output) - 3]
        return output

##############################


class Term(Expression):

    def __init__(self):
        super(Term, self).__init__([self])

    @property
    def terms(self):
        return [self]

    def __mul__(m1, m2):
        assert map(Term.isTerm, (m1, m2)
                   ), "Arguements must be of type ExponentialTerm"
        return MultTerm([m1, m2])

    def simplify(self):
        return self

    @staticmethod
    def isTerm(t):
        return isinstance(t, Term)


class X(Term):

    def __init__(self, var="x"):
        super(X, self).__init__()
        self.var = var

    def evaluate(self, x):
        return x

    def __add__(x1, x2):
        assert map(Expression.isExpression, (x1, x2)
                   ), "Arguements must be an Expression"
        if x1.__class__ == x2.__class__:
            return PowerTerm(2, x1, 1)
        if utils.isnumeric(x2):
            return Expression(x1.terms + [x2])
        return Expression(x1.terms + x2.terms)

    def __mul__(x1, x2):
        assert map(Expression.isExpression, (x1, x2)
                   ), "Arguements must be an Expression"
        if x2.__class__ == Expression:
            return x2.__mul__(x1)
        if x1.__class__ == x2.__class__:
            return PowerTerm(1, x1, 2)
        if x2.__class__ == PowerTerm:
            x2copy = x2.copy()
            return x2copy.set_exp(x2copy.exp + 1)
        if utils.isnumeric(x2):
            return PowerTerm(x2, x1, 1)
        return MultTerm((x1, x2))

    def __neg__(self):
        # ERROR HANDLING DONE IN MUL
        return PowerTerm(-1, self, 1)

    def __str__(self):
        return self.var

    @staticmethod
    def isX(x):
        return isinstance(x, X)

    def derivative(self):
        return 1


class MultTerm(Term):

    def __init__(self, terms):
        super(MultTerm, self).__init__()
        self.termsMultiplied = list(terms)
        self.c = 1
        for term in self.termsMultiplied:
            if utils.isnumeric(term):
                self.c *= term
                self.termsMultiplied.remove(term)

    def __mul__(m1, m2):
        assert map(Expression.isExpression, (m1, m2)
                   ), "Arguements must be of type ExponentialTerm"
        if m2.__class__ == Expression:
            return m2 * m1
        if m1.__class__ == m2.__class__:
            return MultTerm(m1.termsMultiplied + m2.termsMultiplied)
        return MultTerm([m1, m2])

    def simplify(self):
        if self.c == 0:
            return 0
        else:
            self.termsMultiplied = [term.simplify() for term in self.termsMultiplied]

    def evaluate(self, x):
        product = 1
        for term in self.termsMultiplied:
            product *= term.evaluate(x)
        return product

    def __str__(self):
        output = str(self.c) + "*"
        for term in self.termsMultiplied:
            output += "(" + str(term) + ")" + "*"
        return  output[: len(output) - 1]


class PowerTerm(Term):

    def __init__(self, coefficient=1, main=X(), exp=1):
        super(PowerTerm, self).__init__()
        assert utils.isnumeric(coefficient), "Coefficient must be a number"
        assert Expression.isExpression(main)
        assert Expression.isExpression(exp)
        self.c = coefficient
        self.main = main
        self.exp = exp

    def set_exp(self, exp):
        self.exp = exp
        return self

    def simplify(self):
        return self._simplifier()

    def _simplifier(self):
        if self.exp == 0:
            return self.c
        if self.c == 0:
            return 0
        return self

    def __add__(p1, p2):
        assert map(Expression.isExpression, (p1, p2)
                   ), "Arguements must be of type ExponentialTerm"
        if p2 == 0:
            return p1
        if PowerTerm.isPwrTerm(p2) and p1.exp == p2.exp and p1.main == p2.main:
            return PowerTerm(p1.c + p2.c, p1.main, p1.exp)
        else:
            return Expression([p1, p2])

    def __mul__(p1, p2):
        assert map(Expression.isExpression, (p1, p2)
                   ), "Arguements must be of type ExponentialTerm"
        if p2.__class__ == Expression:
            return p2.__mul__(p1)
        if p1.__class__ == p2.__class__ and p1.main == p2.main:
            return PowerTerm(p1.c * p2.c, p1.main, p1.exp + p2.exp)
        if isinstance(p2, X):
            copy = p1.copy()
            return copy.set_exp(copy.exp + 1)
        if utils.isnumeric(p2):
            return PowerTerm(p1.c * p2, p1.main, p1.exp)
        return MultTerm([p1, p2])

    def evaluate(self, x):
        return self.c * pow(self.main(x), self.exp(x))

    @staticmethod
    def isPwrTerm(p):
        return isinstance(p, PowerTerm)

    def __str__(self):
        o_paran = "(" if not Term.isTerm(self.main) else ""
        c_paran = ")" if not Term.isTerm(self.main) else ""

        output = ""

        main = str(self.main)
        exp = "^" + str(self.exp)
        coefficient = str(self.c)

        if self.c == 0:
            return "0"
        if self.exp == 0:
            return str(self.c)
        if self.exp == 1:
            exp = ""
        if self.c == 1:
            return o_paran + main + c_paran + exp
        return o_paran + coefficient + main + c_paran + exp


class ExponentialTerm(Term):

    def __init__(self, coefficient=1, exp=1):
        super(ExponentialTerm, self).__init__()
        assert isinstance(exp, Expression)
        self.c = coefficient
        self.exp = exp

    def __add__(e1, e2):
        assert map(ExponentialTerm.isExpTerm, (e1, e2)
                   ), "Arguements must be of type ExponentialTerm"
        if e2 == 0:
            return e1
        if e1.exp == e2.exp:
            return ExponentialTerm(e1.c + e2.c, e1.exp)
        else:
            return Expression([e1, e2])

    def simplify(self):
        if self.c == 0:
            return 0
        return self

    def __mul__(e1, e2):
        assert map(Expression.isExpression, (e1, e2)
                   ), "Arguements must be of type ExponentialTerm"
        if e2.__class__ == Expression:
            return e2.__mul__(e1)
        if ExponentialTerm.isExpTerm(e2):
            return ExponentialTerm(e1.c * e2.c, e1.exp + e2.exp)
        return MultTerm([e1, e2])

    def evaluate(self, x):
        return self.c * math.exp(self.exp(x))

    def isExpTerm(e):
        return isinstance(e, ExponentialTerm)

    def __str__(self):
        if self.exp == 0:
            return str(self.c)
        if self.c == 0:
            return "0"
        else:
            return str(self.c) + "e" + "^(" + str(self.exp) + "x)"


class LogTerm(Term):
    """c*log_base(Beta * x)"""

    def __init__(self, insideTerm, coefficient=1, base=10):
        super(LogTerm, self).__init__()
        assert Expression.isExpression(insideTerm)
        self.c = coefficient
        self.base = base
        self.insideTerm = insideTerm

    def __mul__(l1, l2):
        assert map(Expression.isExpression, (e1, e2)
                   ), "Arguements must be of type ExponentialTerm"
        if x2.__class__ == Expression:
            return x2.__mul__(x1)
        if LogTerm.isLogTerm(l2) and l2.base == l1.base:
            return LogTerm(PowerTerm(1,l1.insideTerm,l1.c) + PowerTerm(1,l2.insideTerm, l2.c),1, l1.base)
        return MultTerm((l1, l2))

    def simplify(self):
        self.insideTerm = self.insideTerm.simplify()
        return self

    def evaluate(self, x):
        return self.c * math.log(self.insideTerm(x), self.base)

    @staticmethod
    def isLogTerm(l):
        return isinstance(l, LogTerm)

    def __str__(self):
        base = "e" if self.base == math.e else str(self.base)
        coeff = str(self.c) if self.c != 1 else ""
        return coeff + "log_" + base + "(" + str(self.insideTerm) + ")"
