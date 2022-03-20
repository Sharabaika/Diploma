import pandas as pd
import numpy as np
import os

class ResultSaving:
    def __init__(self, *args, **kwargs):
        self.logger = Logger(*args, **kwargs)
        self.params = {}

    def AddParams(self, **kwargs):
        self.params.update(kwargs)

    def SaveResult(self, folder, result_name, category, *args, **kwargs ):
        df = pd.DataFrame(data = kwargs)
        float_format = '{:.5f}'.format

        path = os.path.join(folder, result_name)

        if not os.path.exists(path):
            os.makedirs(path)

        file_path = os.path.join(path, f"{category}.csv")
        df.to_csv(file_path, float_format=float_format)
    
    def SaveResults(self, folder, result_name):
        float_format = '{:.5f}'.format

        path = os.path.join(folder, result_name)

        if not os.path.exists(path):
            os.makedirs(path)

        logger_path = os.path.join(path, "logs.csv")
        self.logger.ToDataFrame().to_csv(logger_path, float_format=float_format)

        if self.params is not None:
            params_path = os.path.join(path, "params.csv")
            pd.DataFrame(data=self.params, index=[0]).to_csv(params_path, float_format=float_format)

    def ReadResult(path):
        df = pd.read_csv(path)
        return df

class Logger:
    def __init__(self, *args, **kwargs):
        self.errors = {}
        for arg in args:
            self.errors[arg] = []

    def LogErrors(self, **kwargs):
        for key, value in kwargs.items():
            self.errors[key].append(value)

    def LogErrorsList(self, **kwargs):
        for key, value in kwargs.items():
            self.errors[key].extend(value)

    def ToDataFrame(self):
        return pd.DataFrame(data=self.errors)


def main():
    logger = Logger("Psi", "W")
    logger.LogErrors(Psi = 1, W = 2)
    logger.LogErrors(Psi = 3, W = 4)
    print(logger.ToDataFrame())


if __name__ == "__main__":
    main()