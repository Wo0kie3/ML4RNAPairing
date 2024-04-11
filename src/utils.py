import numpy as np
import itertools


def _planar_angle(p1, p2, p3):
    b1 = p2 - p1
    b2 = p2 - p3

    angle = np.arccos(
        np.dot(b1, b2) / (np.linalg.norm(b1) * np.linalg.norm(b2)))
    return np.degrees(angle)


def _torsion_angle(p1, p2, p3, p4):
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3

    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)

    torsion = np.arctan2(
        np.dot(np.cross(n1, n2), b2 * np.linalg.norm(b2)), np.dot(n1, n2)
    )
    return np.degrees(torsion)


def distance_3d(p1, p2):
    return np.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2)


def pairwise_distances(coords, num_bins=20):
    distances = np.sqrt(
        np.sum((coords[:, np.newaxis, :] - coords[np.newaxis, :, :]) ** 2, axis=-1))
    upper_triangle_indices = np.triu_indices_from(distances, k=1)
    flattened_distances = distances[upper_triangle_indices]
    histogram, _ = np.histogram(
        flattened_distances, bins=num_bins, density=True)
    normalized_histogram = histogram / np.sum(histogram)
    return normalized_histogram


def p_angles(vector1: list, vector2: list, num_bins: int) -> list:

    pairs_from_vector1 = list(itertools.combinations(vector1, 2))
    pairs_from_vector2 = list(itertools.combinations(vector2, 2))

    # Generate triplets where two elements are from vector1 and one from vector2
    triplets_1 = []
    for single_element in vector2:
        for pair in pairs_from_vector1:
            triplets_1.append((pair[0], single_element, pair[1]))

    # Generate triplets where two elements are from vector2 and one from vector1
    triplets_2 = []
    for single_element in vector1:
        for pair in pairs_from_vector2:
            triplets_2.append((pair[0], single_element, pair[1]))

    # Combine all triplets
    all_triplets = triplets_1 + triplets_2
    angles = []

    all_triplets = np.array(all_triplets)

    angles = [_planar_angle(trip[0], trip[1], trip[2])
              for trip in all_triplets]

    histogram, _ = np.histogram(angles, bins=num_bins, density=True)
    normalized_histogram = histogram / np.sum(histogram)
    return normalized_histogram


def t_angles(vector_i: list, pos_i: int, vector_j: list, pos_j: int, num_bins: int) -> list:

    pairs_from_vector1 = [(vector_i[pos_i], atom)
                          for atom in vector_i if list(atom) != list(vector_i[pos_i])]

    pairs_from_vector2 = [(vector_j[pos_j], atom)
                          for atom in vector_j if list(atom) != list(vector_j[pos_j])]

    combinations_of_pairs = []
    for pair1 in pairs_from_vector1:
        for pair2 in pairs_from_vector2:
            combinations_of_pairs.append(pair1 + pair2)

    combinations_of_pairs = np.array(combinations_of_pairs)

    angles = [_torsion_angle(quad[0], quad[1], quad[2], quad[3])
              for quad in combinations_of_pairs]

    histogram, _ = np.histogram(angles, bins=num_bins, density=True)
    normalized_histogram = histogram / np.sum(histogram)
    return normalized_histogram
