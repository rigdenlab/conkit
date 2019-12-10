from conkit.plot.figure import Figure
import pandas as pd
from bokeh.plotting import figure, show, output_file
from bokeh.models import Range1d, Legend
from bokeh.models.tools import HoverTool
from bokeh.models import ColumnDataSource


class ContactPlotPlusFigure(Figure):

    def __init__(self, seq, conpred, sspred=None,
                 conservationpred=None, mempred = None, dispred=None,
                 data=None, dfObj=None):

        self.conpred = conpred
        self.seq = seq
        self.sspred = sspred
        self.conservationpred = conservationpred
        self.mempred = mempred
        self.dispred = dispred
        self.data = data
        self.dfObj = dfObj


    def contact_data(self):
        # Assign the sequence register to your contact prediction
        self.conpred.sequence = self.seq
        self.conpred.set_sequence_register()

        # tidy contact prediction before plotting
        self.conpred.remove_neighbors(inplace=True)
        self.conpred.sort('raw_score', reverse=True, inplace=True)

        # plot top-L only so slice the contact cmap
        cmap = self.conpred[:self.conpred.sequence.seq_len]

        # place contact and seq data into df
        dfObj = pd.DataFrame(
            columns=['id', 'res1', 'res1_chain', 'res1_seq', 'res2', 'res2_chain', 'res2_seq', 'raw_score'])

        for x in cmap:
            dfObj = dfObj.append({'id': x.id, 'res1': x.res1, 'res1_chain': x.res1_chain, 'res1_seq': x.res1_seq,
                                  'res2': x.res2, 'res2_chain': x.res2_chain, 'res2_seq': x.res2_seq,
                                  'raw_score': x.raw_score}, ignore_index=True)

        return dfObj


    def data_frame_maker(self):
        # create residue list from seq object
        sequence = self.seq
        sequence_for_pred = []
        position = []
        for indx, x in enumerate(sequence.seq):
            sequence_for_pred.append(x)
            position.append(indx + 1)

        # create df to hold predictions and seq & positions
        data = pd.DataFrame({'position': position, 'aa': sequence_for_pred})

        predictions  = {self.sspred:'ss2_prediction', self.conservationpred:'conservation_score',
                        self.mempred:'membrane_prediction', self.dispred:'disorder_prediction'}

        for key, value in predictions.items():
            if key != None:
                temp = key
                temp_pred_list = []
                for y in temp:
                    temp_pred_list.append(getattr(y, value))
                data[value] = temp_pred_list

        return data


    def plot(self):

        #plot contact data
        dfObj = self.contact_data()
        source = ColumnDataSource(dfObj)
        p = figure(x_range=Range1d(0, len(dfObj)), y_range=Range1d(0, len(dfObj)),
                   tools=['pan', 'wheel_zoom', 'box_zoom', 'reset', 'save']
                   , toolbar_location="above",
                   plot_height=575, plot_width=700)
        p1 = p.circle(x='res1_seq', y='res2_seq', size=3, color='black', source=source)
        p.add_tools(HoverTool(renderers=[p1], tooltips=[('Contact', '@id'), ('res1', '@res1'), ('res2', '@res2'),
                                                        ('confidence', '@raw_score')]))
        p1a = p.circle(x='res2_seq', y='res1_seq', size=3, color='black', source=source)
        p.add_tools(HoverTool(renderers=[p1a], tooltips=[('Contact', '@id'), ('res1', '@res1'), ('res2', '@res2'),
                                                         ('confidence', '@raw_score')]))

        #set up legend data
        p2 = p3 = p4 = p5 = p6 = p7 = p8 = p9 = None

        #grap prediction data
        data = self.data_frame_maker()

        #plot mem prediction
        if self.mempred != None:
            data_i = data.loc[(data.membrane_prediction == 'i')]
            xy_i = data_i['position'].tolist()
            data_m = data.loc[(data.membrane_prediction == 'M')]
            xy_m = data_m['position'].tolist()
            data_o = data.loc[(data.membrane_prediction == 'o')]
            xy_o = data_o['position'].tolist()

            x_i = xy_i
            y_i = xy_i
            x_m = xy_m
            y_m = xy_m
            x_o = xy_o
            y_o = xy_o

            p2 = p.circle(x_i, y_i, size=10, color='green')
            p3 = p.circle(x_m, y_m, color='red', size=10)
            p4 = p.circle(x_o, y_o, color='yellow', size=10)

        #plot ss prediction
        if self.sspred != None:
            #create a seperate data set for each ss classification
            data_h = data.loc[(data.ss2_prediction == 'H')]
            xy_h = data_h['position'].tolist()
            xy_h_aa_name = data_h['aa'].tolist()
            data_c = data.loc[(data.ss2_prediction == 'C')]
            xy_c = data_c['position'].tolist()
            xy_c_aa_name = data_c['aa'].tolist()
            data_e = data.loc[(data.ss2_prediction == 'E')]
            xy_e = data_e['position'].tolist()
            xy_e_aa_name = data_e['aa'].tolist()

            # initialise aa property dictionary
            aa_properites = {'A': 'Non-polar, aliphatic',
                             'R': 'Positively charged (basic; non-acidic); Polar; Hydrophilic; pK=12.5',
                             'N': 'Polar, non-charged',
                             'D': 'Negatively charged (acidic); Polar; Hydrophilic; pK=3.9',
                             'C': 'Polar, non-charged',
                             'E': 'Negatively charged (acidic); Polar; Hydrophilic; pK=4.2',
                             'Q': 'Polar, non-charged',
                             'G': 'Non-polar, aliphatic',
                             'H': 'Positively charged (basic; non-acidic); Polar; Hydrophilic; pK=6.0',
                             'I': 'Non-polar, aliphatic',
                             'L': 'Non-polar, aliphatic',
                             'K': 'Positively charged (basic; non-acidic); Polar; Hydrophilic; pK=10.5',
                             'M': 'Polar, non-charged',
                             'F': 'Aromatic',
                             'P': 'Non-polar, aliphatic',
                             'S': 'Polar, non-charged',
                             'T': 'Polar, non-charged',
                             'W': 'Aromatic',
                             'Y': 'Aromatic',
                             'V': 'Aromatic'
                             }

            #add aa properties to each residue from aa property dictionary
            cnt = 0
            h_aa_properies = []
            for x in xy_h_aa_name:
                for y in aa_properites:
                    if x == y:
                        h_aa_properies.append(aa_properites[y])
                        cnt += 1

            c_aa_properies = []
            for x in xy_c_aa_name:
                for y in aa_properites:
                    if x == y:
                        c_aa_properies.append(aa_properites[y])

            e_aa_properies = []
            for x in xy_e_aa_name:
                for y in aa_properites:
                    if x == y:
                        e_aa_properies.append(aa_properites[y])


            #plot ss on to the contact plot with residue properties
            source2 = ColumnDataSource(data=dict(x=xy_h, y=xy_h, aa_name=xy_h_aa_name, aa_properties=h_aa_properies))
            source3 = ColumnDataSource(data=dict(x=xy_c, y=xy_c, aa_name=xy_c_aa_name, aa_properties=c_aa_properies))
            source4 = ColumnDataSource(data=dict(x=xy_e, y=xy_e, aa_name=xy_e_aa_name, aa_properties=e_aa_properies))

            p5 = p.circle(x='x', y='y', size=3, color='orange', source=source2)
            p.add_tools(HoverTool(renderers=[p5], tooltips=[('Position', '@x'), ('Residue Name', '@aa_name'),
                                                            ('Properties', '@aa_properties')]))
            p6 = p.circle(x='x', y='y', color='blue', size=3, source=source3)
            p.add_tools(HoverTool(renderers=[p6], tooltips=[('Position', '@x'), ('Residue Name', '@aa_name'),
                                                            ('Properties', '@aa_properties')]))
            p7 = p.circle(x='x', y='y', color='pink', size=3, source=source4)
            p.add_tools(HoverTool(renderers=[p7], tooltips=[('Position', '@x'), ('Residue Name', '@aa_name'),
                                                            ('Properties', '@aa_properties')]))

        # plot consurf values
        if self.conservationpred != None:
            consurf_values = list(range(1, 10))

            x_dct = {}
            for x in consurf_values:
                conserved_aa = data.loc[(data.conservation_score == x)]
                x_dct[x] = conserved_aa['position'].tolist()

            y_dct = {}
            for y in consurf_values:
                conserved_aa = data.loc[(data.conservation_score == y)]
                y_dct[y] = conserved_aa['position'].tolist()
                y_dct[y] = [int(x) + 3 for x in y_dct[y]]

            Blues = {1: '#f7fbff', 2: '#deebf7', 3: '#c6dbef', 4: '#9ecae1', 5: '#6baed6', 6: '#4292c6',
                         7: '#2171b5', 8: '#08519c', 9: '#08306b'}

            for x in Blues:
                p.triangle(x_dct[x], y_dct[x], size=5, color=Blues[x])

        # plot disorder prediction
        if self.dispred != None:
            data_disorder = data.loc[(data.disorder_prediction >= 0.5)]
            x_d = data_disorder['position'].tolist()
            y_d = []
            for x in x_d:
                x = int(x) - 3
                y_d.append(x)

            data_not_disorder = data.loc[(data.disorder_prediction < 0.5)]
            x_not_dis = data_not_disorder['position'].tolist()
            y_not_dis = []
            for x in x_not_dis:
                x = int(x) - 3
                y_not_dis.append(x)

            x_d = x_d
            y_d = y_d
            x_nd = x_not_dis
            y_nd = y_not_dis

            p8 = p.circle(x_d, y_d, size=5, color='lime')
            p9 = p.circle(x_nd, y_nd, color='grey', size=5)


        #construct legend
        legend_data = {"Mem_pred-I": [p2], "Mem_pred-M": [p3], "Mem_pred-O": [p4],
                      "SS-H": [p5], "SS-C": [p6], "SS-E": [p7], "Disorder": [p8], "Ordered": [p9]}

        items = [("Contact", [p1])]

        for key, value in legend_data.items():
            if value != [None]:
                legend_item = (key,value)
                items.append(legend_item)

        legend = Legend(items=items, location="center")

        p.add_layout(legend, 'right')

        #save and show plot
        output_file('/Users/shahrammesdaghi/Downloads/test.html')
        show(p)
        #return output


if __name__ == '__main__':

    import conkit.io

    conpred = conkit.io.read('/Users/shahrammesdaghi/Downloads/w9dy28_meta_respre.evfold', 'evfold').top
    seq = conkit.io.read('/Users/shahrammesdaghi/Downloads/w9dy28.fasta', 'fasta').top
    sspred = conkit.io.read('/Users/shahrammesdaghi/Downloads/w9dy28.ss2', 'psipred')
    conservationpred = conkit.io.read('/Users/shahrammesdaghi/Downloads/consurf.grades.txt', 'consurf')
    mempred = conkit.io.read('/Users/shahrammesdaghi/Downloads/query.result.txt', 'topcons')
    dispred = conkit.io.read('/Users/shahrammesdaghi/Downloads/iupred2a.txt', 'iupred2a')

    #test_plot = ContactPlotPlusFigure(seq, conpred, ss_prediction, con_pred, mem_pred, dis_pred)
    #test_plot.plot_plus()

    test_plot_2 = conkit.plot.ContactPlotPlusFigure(seq=seq, conpred=conpred, conservationpred=conservationpred,
                                                    mempred=mempred, dispred=dispred, sspred=sspred)
    test_plot_2.plot()


