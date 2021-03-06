from anastruct.fem.elements import det_moment, det_shear
import numpy as np
import math

from typing import TYPE_CHECKING, List, Tuple

if TYPE_CHECKING:
    from anastruct.fem.system import SystemElements
    from anastruct.fem.elements import Element


def set_force_vector(
    system: "SystemElements", force_list: List[Tuple[int, int, float]]
):
    """
    :param force_list: list containing tuples with the
    1. number of the node,
    2. the number of the direction (1 = x, 2 = z, 3 = y)
    3. the force
    [(1, 3, 1000)] node=1, direction=3 (y), force=1000
    list may contain multiple tuples
    :return: Vector with forces on the nodes
    """
    assert system.system_force_vector is not None
    for id_, direction, force in force_list:
        system.system_force_vector[(id_ - 1) * 3 + direction - 1] += force

    return system.system_force_vector


def prep_matrix_forces(system: "SystemElements"):
    system.system_force_vector = system.system_force_vector = np.zeros(
        len(system._vertices) * 3
    )
    apply_perpendicular_q_load(system)
    apply_point_load(system)
    apply_moment_load(system)


def apply_moment_load(system: "SystemElements"):
    for node_id, Ty in system.loads_moment.items():
        set_force_vector(system, [(node_id, 3, Ty)])


def apply_point_load(system: "SystemElements"):
    for node_id in system.loads_point:
        Fx, Fz = system.loads_point[node_id]
        # system force vector.
        set_force_vector(
            system,
            [
                (node_id, 1, Fx),
                (node_id, 2, Fz),
            ],
        )


def apply_perpendicular_q_load(system: "SystemElements"):
    for element_id in system.loads_dead_load:
        element = system.element_map[element_id]
        if element.q_load is None and element.dead_load == 0:
            continue

        qi_perpendicular = element.all_q_load[0]
        q_perpendicular = element.all_q_load[1]
        if not (
            math.isclose(element.q_load[0] + element.dead_load, qi_perpendicular)
        ) and not (
            math.isclose(element.q_load[1] + element.dead_load, q_perpendicular)
        ):
            apply_parallel_q_load(system, element)
        if q_perpendicular == 0 and qi_perpendicular == 0:
            continue

        kl = element.constitutive_matrix[1][1] * 1e6
        kr = element.constitutive_matrix[2][2] * 1e6

        # minus because of systems positive rotation
        left_moment = det_moment(
            kl, kr, qi_perpendicular, q_perpendicular, 0, element.EI, element.l
        )
        right_moment = -det_moment(
            kl, kr, qi_perpendicular, q_perpendicular, element.l, element.EI, element.l
        )
        rleft = det_shear(
            kl, kr, qi_perpendicular, q_perpendicular, 0, element.EI, element.l
        )
        rright = -det_shear(
            kl, kr, qi_perpendicular, q_perpendicular, element.l, element.EI, element.l
        )

        rleft_x = rleft * math.sin(element.a1)
        rright_x = rright * math.sin(element.a2)

        rleft_z = rleft * math.cos(element.a1)
        rright_z = rright * math.cos(element.a2)

        if element.type == "truss":
            left_moment = 0
            right_moment = 0

        primary_force = np.array(
            [rleft_x, rleft_z, left_moment, rright_x, rright_z, right_moment]
        )
        element.element_primary_force_vector -= primary_force

        # Set force vector
        assert system.system_force_vector is not None
        system.system_force_vector[
            (element.node_id1 - 1) * 3 : (element.node_id1 - 1) * 3 + 3
        ] += primary_force[0:3]
        system.system_force_vector[
            (element.node_id2 - 1) * 3 : (element.node_id2 - 1) * 3 + 3
        ] += primary_force[3:]


def apply_parallel_q_load(system: "SystemElements", element: "Element"):
    direction = element.q_direction
    # dead load
    factor_dl = abs(math.sin(element.angle))

    def update(Fx_left: float, Fx_right, Fz_left: float, Fz_right):
        element.element_primary_force_vector[0] -= Fx_left
        element.element_primary_force_vector[1] -= Fz_left
        element.element_primary_force_vector[3] -= Fx_right
        element.element_primary_force_vector[4] -= Fz_right

        set_force_vector(
            system,
            [
                (element.node_1.id, 2, Fz_left),
                (element.node_2.id, 2, Fz_right),
                (element.node_1.id, 1, Fx_left),
                (element.node_2.id, 1, Fx_right),
            ],
        )

    if direction == "x":
        factor = abs(math.cos(element.angle))

        q_load = [i * factor for i in element.q_load]
        for q_element in (q_load, element.dead_load * factor_dl):
            # q_load working at parallel to the elements x-axis          # set the proper direction
            if isinstance(q_element, list):
                qi = q_element[0]
                q = q_element[1]
            else:
                qi = q = q_element

            if q == 0 and qi == 0:
                cg = 0.5
            else:
                cg = (element.l / 3) * ((qi + 2 * q) / (qi + q))

            eq = ((q + qi) * element.l) / 2

            Fx_left = -eq * math.cos(element.angle) * (element.l - cg) / element.l

            Fx_right = -eq * math.cos(element.angle) * cg / element.l

            Fz_left = (
                eq
                * abs(math.sin(element.angle))
                * (element.l - cg)
                / element.l
                * np.sign(math.sin(element.angle))
            )

            Fz_right = (
                eq
                * abs(math.sin(element.angle))
                * cg
                / element.l
                * np.sign(math.sin(element.angle))
            )

            update(Fx_left, Fx_right, Fz_left, Fz_right)

    else:
        if math.isclose(element.angle, 0):
            # horizontal element cannot have parallel forces due to self weight or q-load in y direction.
            return None
        factor = abs(math.sin(element.angle))

        q_load = [i * factor for i in element.q_load]
        for q_element in (q_load, element.dead_load * factor_dl):
            # q_load working at parallel to the elements x-axis          # set the proper direction
            if isinstance(q_element, list):
                qi = q_element[0]
                q = q_element[1]
            else:
                qi = q = q_element

            if q == 0 and qi == 0:
                cg = 0.5
            else:
                cg = (element.l / 3) * ((qi + 2 * q) / (qi + q))

            eq = ((q + qi) * element.l) / 2

            Fx_left = (
                eq
                * math.cos(element.angle)
                * (element.l - cg)
                / element.l
                * -np.sign(math.sin(element.angle))
            )

            Fx_right = (
                eq
                * math.cos(element.angle)
                * cg
                / element.l
                * -np.sign(math.sin(element.angle))
            )

            Fz_left = eq * abs(math.sin(element.angle)) * (element.l - cg) / element.l

            Fz_right = eq * abs(math.sin(element.angle)) * cg / element.l

            update(Fx_left, Fx_right, Fz_left, Fz_right)


def dead_load(system: "SystemElements", g: float, element_id: int):
    system.loads_dead_load.add(element_id)
    system.element_map[element_id].dead_load = g


def assemble_system_matrix(
    system: "SystemElements", validate: bool = False, geometric_matrix: bool = False
):
    """
    Shape of the matrix = n nodes * n d.o.f.
    Shape = n * 3
    """
    system._remainder_indexes = []
    if not geometric_matrix:
        shape = len(system.node_map) * 3
        system.shape_system_matrix = shape
        system.system_matrix = np.zeros((shape, shape))

    assert system.system_matrix is not None
    for matrix_index, K in system.system_spring_map.items():
        #  first index is row, second is column
        system.system_matrix[matrix_index][matrix_index] += K

    # Determine the elements location in the stiffness matrix.
    # system matrix [K]
    #
    # [fx 1] [K         |  \ node 1 starts at row 1
    # |fz 1] |  K       |  /
    # |Ty 1] |   K      | /
    # |fx 2] |    K     |  \ node 2 starts at row 4
    # |fz 2] |     K    |  /
    # |Ty 2] |      K   | /
    # |fx 3] |       K  |  \ node 3 starts at row 7
    # |fz 3] |        K |  /
    # [Ty 3] [         K] /
    #
    #         n   n  n
    #         o   o  o
    #         d   d  d
    #         e   e  e
    #         1   2  3
    #
    # thus with appending numbers in the system matrix: column = row

    for i in range(len(system.element_map)):
        element = system.element_map[i + 1]
        element_matrix = element.stiffness_matrix

        # n1 and n2 are starting indexes of the rows and the columns for node 1 and node 2
        n1 = (element.node_1.id - 1) * 3
        n2 = (element.node_2.id - 1) * 3
        system.system_matrix[n1 : n1 + 3, n1 : n1 + 3] += element_matrix[0:3, :3]
        system.system_matrix[n1 : n1 + 3, n2 : n2 + 3] += element_matrix[0:3, 3:]

        system.system_matrix[n2 : n2 + 3, n1 : n1 + 3] += element_matrix[3:6, :3]
        system.system_matrix[n2 : n2 + 3, n2 : n2 + 3] += element_matrix[3:6, 3:]

    # returns True if symmetrical.
    if validate:
        assert np.allclose((system.system_matrix.transpose()), system.system_matrix)


def set_displacement_vector(system, nodes_list):
    """
    :param nodes_list: list containing tuples with
    1.the node
    2. the d.o.f. that is set
    :return: Vector with the displacements of the nodes (If displacement is not known, the value is set
    to NaN)
    """
    if system.system_displacement_vector is None:
        system.system_displacement_vector = np.ones(len(system._vertices) * 3) * np.NaN

    for i in nodes_list:
        index = (i[0] - 1) * 3 + i[1] - 1

        try:
            system.system_displacement_vector[index] = 0
        except IndexError as e:
            raise IndexError(
                e,
                "This often occurs if you set supports before the all the elements are modelled. "
                "First finish the model.",
            )
    return system.system_displacement_vector


def process_conditions(system):
    indexes = []
    # remove the unsolvable values from the matrix and vectors
    for i in range(system.shape_system_matrix):
        if system.system_displacement_vector[i] == 0:
            indexes.append(i)
        else:
            system._remainder_indexes.append(i)

    system.system_displacement_vector = np.delete(
        system.system_displacement_vector, indexes, 0
    )
    system.reduced_force_vector = np.delete(system.system_force_vector, indexes, 0)
    system.reduced_system_matrix = np.delete(system.system_matrix, indexes, 0)
    system.reduced_system_matrix = np.delete(system.reduced_system_matrix, indexes, 1)


def process_supports(system):
    for node in system.supports_hinged:
        set_displacement_vector(system, [(node.id, 1), (node.id, 2)])

    for i in range(len(system.supports_roll)):
        if not system.supports_roll_rotate[i]:
            set_displacement_vector(
                system,
                [
                    (system.supports_roll[i].id, system.supports_roll_direction[i]),
                    (system.supports_roll[i].id, 3),
                ],
            )
        else:
            set_displacement_vector(
                system,
                [(system.supports_roll[i].id, system.supports_roll_direction[i])],
            )

    for node in system.supports_rotational:
        set_displacement_vector(system, [(node.id, 3)])

    for node in system.internal_hinges:
        set_displacement_vector(system, [(node.id, 3)])
        for el in system.node_element_map[node.id]:
            el.compile_constitutive_matrix(
                el.EA, el.EI, el.l, el.springs, el.node_1.hinge, el.node_2.hinge
            )
            el.compile_stiffness_matrix()

    for node in system.supports_fixed:
        set_displacement_vector(system, [(node.id, 1), (node.id, 2), (node.id, 3)])

    for node, roll in system.supports_spring_x:
        if not roll:
            set_displacement_vector(system, [(node.id, 2)])

    for node, roll in system.supports_spring_z:
        if not roll:
            set_displacement_vector(system, [(node.id, 1)])

    for node, roll in system.supports_spring_y:
        if not roll:
            set_displacement_vector(system, [(node.id, 1), (node.id, 2)])

    for node_id, angle in system.inclined_roll.items():
        for el in system.node_element_map[node_id]:
            if el.node_1.id == node_id:
                el.a1 = el.angle + angle
            elif el.node_2.id == node_id:
                el.a2 = el.angle + angle

            el.compile_kinematic_matrix(el.a1, el.a2, el.l)
            el.compile_stiffness_matrix()
