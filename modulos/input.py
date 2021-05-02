# Module to read and organize the input Excel spreadsheet.


class inputInfo():

    # filename = string containing the excel input filename location
    def __init__(self, filename):
        import pandas as pd

        try:
            headerList = pd.read_excel(
                filename,
                sheet_name='Input',
                skiprows=1,
                usecols='P',
                nrows=3,
                header=None,
                # dtype=object
                ).T

            projectInfo = pd.read_excel(
                filename,
                sheet_name='Input',
                skiprows=5,
                usecols='B,E,G,H,J,L,O',
                nrows=1,
                # dtype=object
                )

            filledData = pd.read_excel(
                filename,
                sheet_name='Input',
                skiprows=10,
                usecols='B:P',
                # dtype=object
                )

            heatFlux = pd.read_excel(
                filename,
                sheet_name='Coefs',
                skiprows=2,
                usecols='C',
                # dtype=object
                )

            convCoef = pd.read_excel(
                filename,
                sheet_name='Coefs',
                skiprows=2,
                usecols='D',
                # dtype=object
                )

        # except FileNotFoundError:
        #     raise InputReadError(filename)

        except:
            # Exception to get possible modules not found
            raise

        self.header = headerList.values.tolist()[0]
        self.project = projectInfo.values.tolist()[0]
        self.read = filledData.values.tolist()[0]
        self.q = heatFlux.values.tolist()
        self.h = convCoef.values.tolist()


    # Retorna os parametros da simulacao
    def getParam(self):
        return self.read[1:5]

    # Retorna os inputs do disco
    def getGeom(self):
        return self.read[5:10]

    # Retorna as propriedades do disco
    def getProperties(self):
        return self.read[10:13]

    # Retorna as temperaturas inicial e ambiente
    def getAmb(self):
        return self.read[13:15]

    # Retorna as temperaturas inicial e ambiente
    def getCoefs(self,coef):
        if coef == 'q':
            return self.q
        elif coef == 'h':
            return self.h

    # Retorna as informacoes do projeto
    def projectInfo(self):
        return self.project

    # Retorna o numero de ciclos
    # Para a frenagem unica = 1
    # Para multiplas frenagens > 1
    def ciclos(self):
        if self.project[3] == "Unica":
            return 1
        else:
            return self.project[4]