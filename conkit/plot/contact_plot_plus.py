import conkit.io
import conkit.plot
from conkit.plot.figure import Figure
import pandas as pd
from bokeh.plotting import figure, show, output_file
from bokeh.models import Range1d, Legend
from bokeh.models.tools import HoverTool
from bokeh.models import ColumnDataSource

class ContactPlotPlus(Figure):

    def __init__(self, seq, conpred, ss_prediction = None, con_pred = None, mem_pred = None, dis_pred = None):

        self.conpred = conpred
        self.seq = seq
        self.ss_prediction = ss_prediction
        self.con_pred = con_pred
        self.mem_pred = mem_pred
        self.dis_pred = dis_pred

    def plot_plus(self):

        # Assign the sequence register to your contact prediction
        self.conpred.sequence = self.seq
        self.conpred.set_sequence_register()

        # tidy contact prediction before plotting
        self.conpred.remove_neighbors(inplace=True)
        self.conpred.sort('raw_score', reverse=True, inplace=True)

        # plot top-L only so slice the contact cmap
        cmap = self.conpred[:self.conpred.sequence.seq_len]

        #place contact and seq data into df
        dfObj = pd.DataFrame(
            columns=['id', 'res1', 'res1_chain', 'res1_seq', 'res2', 'res2_chain', 'res2_seq', 'raw_score'])

        for x in cmap:
            dfObj = dfObj.append({'id': x.id, 'res1': x.res1, 'res1_chain': x.res1_chain, 'res1_seq': x.res1_seq,
                                  'res2': x.res2, 'res2_chain': x.res2_chain, 'res2_seq': x.res2_seq,
                                  'raw_score': x.raw_score}, ignore_index=True)

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
        #initialise aa property dictionary
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

        # create residue list from seq object
        sequence = self.seq
        sequence_for_pred = []
        position =[]
        for indx, x in enumerate(sequence.seq):
            sequence_for_pred.append(x)
            position.append(indx+1)

        #create df to hold predictions and seq
        data = pd.DataFrame({'position': position,
                            'aa': sequence_for_pred})

        # plot membrane prediction data if present
        if self.mem_pred != None:
            mem_pred = self.mem_pred
            membrane_pred_list = []
            #position = []
            for x in mem_pred:
                membrane_pred_list.append(x.membrane_prediction)
                #position.append(x.res_seq)

            #data['position'] = position
            data['mem_pred'] = membrane_pred_list
            zeros = [0] * len(membrane_pred_list)
            data['zeros'] = zeros
            data_i = data.loc[(data.mem_pred == 'i')]
            xy_i = data_i['position'].tolist()
            data_m = data.loc[(data.mem_pred == 'M')]
            xy_m = data_m['position'].tolist()
            data_o = data.loc[(data.mem_pred == 'o')]
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

        #hold ss pred in df
        if self.ss_prediction != None:
            ss_pred = self.ss_prediction
            ss_pred_list = []
            for x in ss_pred:
                ss_pred_list.append(x.ss2_prediction)

            data['ss_prediction'] = ss_pred_list

            #create a seperate data set for each ss prediction
            data_h = data.loc[(data.ss_prediction == 'H')]
            xy_h = data_h['position'].tolist()
            xy_h_aa_name = data_h['aa'].tolist()
            data_c = data.loc[(data.ss_prediction == 'C')]
            xy_c = data_c['position'].tolist()
            xy_c_aa_name = data_c['aa'].tolist()
            data_e = data.loc[(data.ss_prediction == 'E')]
            xy_e = data_e['position'].tolist()
            xy_e_aa_name = data_e['aa'].tolist()

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

        #plot consurf data if present
        # TODO is con_pred correct??????
        if self.con_pred != None:
            con_pred = self.con_pred
            cons_score = []
            for x in con_pred:
                cons_score.append(x.conservation_score)

            data['consurf_pred'] = cons_score

            consurf_values = list(range(1, 10))
            x_dct = {}

            for x in consurf_values:
                conserved_aa = data.loc[(data.consurf_pred == x)]
                x_dct[x] = conserved_aa['position'].tolist()

            y_dct = {}

            for y in consurf_values:
                conserved_aa = data.loc[(data.consurf_pred == y)]
                y_dct[y] = conserved_aa['position'].tolist()
                y_dct[y] = [int(x) + 3 for x in y_dct[y]]

            Blues = {1: '#f7fbff', 2: '#deebf7', 3: '#c6dbef', 4: '#9ecae1', 5: '#6baed6', 6: '#4292c6',
                         7: '#2171b5', 8: '#08519c', 9: '#08306b'}

            for x in Blues:
                p.triangle(x_dct[x], y_dct[x], size=5, color=Blues[x])

        # plot disorder prediction data if present
        if self.dis_pred != None:
            dis_pred = self.dis_pred
            dis_pred_list = []
            for x in dis_pred:
                dis_pred_list.append(x.disorder_prediction)

            data['disorder_pred'] = dis_pred_list

            data_disorder = data.loc[(data.disorder_pred >= 0.5)]
            x_d = data_disorder['position'].tolist()
            y_d = []
            for x in x_d:
                x = int(
                    x) - 3  #############################adjust depending on size of prot eg 168 aa -3, 404 aa =7, 236=5
                y_d.append(x)

            data_not_disorder = data.loc[(data.disorder_pred < 0.5)]
            x_not_dis = data_not_disorder['position'].tolist()
            y_not_dis = []
            for x in x_not_dis:
                x = int(
                    x) - 3  #############################adjust depending on size of prot eg 168 aa -3, 404 aa =7, 236=5
                y_not_dis.append(x)

            x_d = x_d
            y_d = y_d
            x_nd = x_not_dis
            y_nd = y_not_dis

            p8 = p.circle(x_d, y_d, size=5, color='lime')
            p9 = p.circle(x_nd, y_nd, color='grey', size=5)

        legend = Legend(items=[
            ("Contact", [p1]),
            ("Mem_pred-I", [p2]), ("Mem_pred-M", [p3]), ("Mem_pred-O", [p4]),
            ("SS-H", [p5]), ("SS-C", [p6]), ("SS-E", [p7]),
            ("Disorder", [p8]), ("Ordered", [p9])
        ], location="center")

        p.add_layout(legend, 'right')
        show(p)


        #save and show plot
        output = output_file('/Users/shahrammesdaghi/Downloads/test.html')
        show(p)
        #return output



if __name__ == '__main__':
    conpred = conkit.io.read('/Users/shahrammesdaghi/Downloads/w9dy28_meta_respre.evfold', 'evfold').top
    seq = conkit.io.read('/Users/shahrammesdaghi/Downloads/w9dy28.fasta', 'fasta').top
    ss_prediction = conkit.io.read('/Users/shahrammesdaghi/Downloads/w9dy28.ss2', 'psipred')
    con_pred = conkit.io.read('/Users/shahrammesdaghi/Downloads/consurf.grades', 'consurf')
    mem_pred = conkit.io.read('/Users/shahrammesdaghi/Downloads/query.result.txt', 'topcons')
    dis_pred = conkit.io.read('/Users/shahrammesdaghi/Downloads/iupred2a.txt', 'iupred2a')

    test_plot = ContactPlotPlus(seq, conpred, ss_prediction, con_pred, mem_pred, dis_pred)
    test_plot.plot_plus()

    #conkit.plot.ContactPlotPlusFigure(conpred, seq, ss_prediction, con_pred, mem_pred, dis_pred)

