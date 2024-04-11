import pandas as pd
import pickle

from utils import gen_quad_and_angle, gen_triplets_and_angle, pairwise_distances


class RNAParing:
    def __init__(self) -> None:
        self.pair_classifier = None
        self.stacking_classifier = None

    def read_csv(self, path):
        """
        Reads a CSV file containing RNA structure data and returns it as a DataFrame.
        """
        return pd.read_csv(path)

    def find_stackings(self, data):
        """
        Processes the data to find stackings. 
        Returns a list of data suitable for classification.
        """
        return []

    def find_pairs(self, data):
        """
        Processes the data to find pairs.
        Returns a list of data suitable for classification.
        """
        return []

    def data_to_representation(self, data):
        pairwise, _ = pairwise_distances()
        planar_angles, _ = gen_triplets_and_angle()
        torsion_angles, _ = gen_quad_and_angle()
        return [*list(pairwise), *list(torsion_angles), *list(planar_angles)]

    def prepare_data_for_classification(self, data):
        """
        Uses find_stackings and find_pairs to prepare data for the classifiers.
        Returns two lists, each containing 150 values derived from the prepared data.
        """

        matrix_nested_dict_names = data.groupby(['file', 'nt']).progress_apply(
            lambda group: group.apply(
                lambda row: [row['atom'], row['x'], row['y'], row['z']], axis=1).tolist()
        ).to_dict()

        stackings = self.find_stackings(data)
        pairs = self.find_pairs(data)

        
        stacking_data = self.data_to_representation(data)  # Derived from stackings/pairs
        pair_data = self.data_to_representation(data)  # Derived from stackings/pairs

        return stacking_data, pair_data
    
    def predict(self, filepath):
        pass
