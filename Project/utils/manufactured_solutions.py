import re

analytical_solutions = {
    'u': '1/2 * r * sin(theta)',
    'v': 'r * cos(theta)',
    'w': '3/2 * cos(theta)',
    'lapl_x': '3/2 * sin(theta) / r',
    'lapl_y': '0*r',
    'h_x': 'r/4 * sin(theta)**2 - r/2 * cos(theta)**2',
    'h_y': '0*r',
    'p': '-1/8 * r * r * sin(theta)**2',
    'grad_p_x': '-1/4 * r * sin(theta)**2',
    'grad_p_y': '-1/4 * r * sin(theta) * cos(theta)',
    'u_star': '{u} + dt * (-{h_x} - {grad_p_x} + nu * {lapl_x})',
    'v_star': '{v} + dt * (-{h_y} - {grad_p_y} + nu * {lapl_y})',
    'rhs' : '1/4 * cos(2 * theta) + 1'
}

pattern = r'{((?:[^\\{}]+|\\.)*)}'

def repl(matchobj):
    return analytical_solutions[matchobj.group(1)]

def parse_solution(expression):
    if re.search(pattern, expression):
        expression = re.sub(pattern, repl, expression)

        if re.search(pattern, expression):
            return parse_solution(expression)
        else:
            return expression
    else:
        return expression