import re

"""
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
}"""

analytical_solutions = {
    'u': '(1 - R**2/ (r**2)) * cos(theta)',
    'v': '- (1 + R**2 / (r**2)) * sin(theta)',
    'w': '0*r',
    'lapl_x': '0*r',
    'lapl_y': '0*r',
    'h_x': '2 * R**2 * cos(2 * theta) / (r**3) - 2 * R**4 / r**5',
    'h_y': '2 * R**2 * sin(2 * theta) / (r**3)',
    'p': '-1/2 * (1 - R**2 / r**2)**2 * cos(theta)**2',
    'grad_p_x': '2 * R**2 * cos(theta) **2 * (R**2 - r**2) / r**5',
    'grad_p_y': '-sin(2 * theta) * (r**4 - 2 * R**2 * r**2 + 2 * R**4) / r**5',
    'u_star': '{u} + dt * (-{h_x} - {grad_p_x} + nu * {lapl_x})',
    'v_star': '{v} + dt * (-{h_y} - {grad_p_y} + nu * {lapl_y})',
    'rhs' : '- 8 * R**4 / r**6 + 2 * R**2 * cos(theta)**2 / r + (R**2 - r**2) * (2 * cos(theta)**2 - 1) / r**6'
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