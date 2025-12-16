# START
p1, p2, q = var("p1 p2 q")
a, b = var("a b")
w = var("w", latex_name=r"\omega")
assume(p1 > 0, p2 > 0, q > 0, a < b)

ga = lambda p: exp(-p ^ 2)
p3 = p1 - w
p4 = p2 + w
wMin = max_symbolic(-q, q - 2 * p2)
wMax = min_symbolic(q, 2 * p1 - q)
Statistical = ga(p2)


integral_result = integrate(Statistical / (2 * pi), w, a, b).simplify_full()
print("After w integration:", integral_result)

# p1 >= p2 case
# q integral split into 3 parts
integralq_1 = integrate(integral_result.subs({a:-q         , b: q         }), q, 0, p2).simplify_full()
integralq_2 = integrate(integral_result.subs({a: q - 2 * p2, b: q         }), q, p2, p1).simplify_full()
integralq_3 = integrate(integral_result.subs({a: q - 2 * p2, b: 2 * p1 - q}), q, p1, p1 + p2).simplify_full()
integral_result1 = (integralq_1 + integralq_2 + integralq_3).simplify_full()
# p1 < p2 case
# q integral split into 3 parts
integralq_1 = integrate(integral_result.subs({a:-q         , b: q         }), q, 0, p1).simplify_full()
integralq_2 = integrate(integral_result.subs({a:-q         , b: 2 * p1 - q}), q, p1, p2).simplify_full()
integralq_3 = integrate(integral_result.subs({a: q - 2 * p2, b: 2 * p1 - q}), q, p2, p1 + p2).simplify_full()
integral_result2 = (integralq_1 + integralq_2 + integralq_3).simplify_full()
print("After q integration:", integral_result1, " for p1>=p2, ", integral_result2, " for p1<p2")
integral_result1 = integrate(integral_result1, p2, 0, p1).simplify_full()
integral_result2 = integrate(integral_result2, p2, p1, oo).simplify_full()
integral_result = (integral_result1 + integral_result2).simplify_full()
print("After p2 integration:", integral_result)

# END


# Copy result to clipboard
# check if xclip is installed
import os

if os.system("which xclip > /dev/null 2>&1") == 0:
    os.system(
        "echo '{}' | xclip -selection clipboard".format(latex(integral_result))
    )
