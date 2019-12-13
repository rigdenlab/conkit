from conkit.plot.figure import Figure
from bokeh.plotting import figure, show, output_file
from bokeh.models import Range1d, Legend
from bokeh.models.tools import HoverTool
from bokeh.models import ColumnDataSource


class ContactPlotPlusFigure(Figure):

    def __init__(self, seq, conpred, sspred=None,
                 conservationpred=None, mempred = None, dispred=None):

        self.conpred = conpred
        self.seq = seq
        self.sspred = sspred
        self.conservationpred = conservationpred
        self.mempred = mempred
        self.dispred = dispred


    def plot(self):
        #plot contact data
        id = []
        res1 = []
        res1_chain = []
        res1_seq = []
        res2 = []
        res2_chain = []
        res2_seq = []
        raw_score = []

        # Assign the sequence register to your contact prediction
        self.conpred.sequence = self.seq
        self.conpred.set_sequence_register()

        # tidy contact prediction before plotting
        self.conpred.remove_neighbors(inplace=True)
        self.conpred.sort('raw_score', reverse=True, inplace=True)

        # plot top-L only so slice the contact cmap
        cmap = self.conpred[:self.conpred.sequence.seq_len]

        for x in cmap:
            id.append(x.id)
            res1.append(x.res1)
            res1_chain.append(x.res1_chain)
            res1_seq.append(x.res1_seq)
            res2.append(x.res2)
            res2_chain.append(x.res2_chain)
            res2_seq.append(x.res2_seq)
            raw_score.append(x.raw_score)

        source = ColumnDataSource(data=dict(id=id, res1=res1, res1_chain=res1_chain,res1_seq=res1_seq,
                                            res2=res2, res2_chain=res2_chain,res2_seq=res2_seq, raw_score=raw_score))

        p = figure(x_range=Range1d(0, len(cmap)), y_range=Range1d(0, len(cmap)),
                   tools=['pan', 'wheel_zoom', 'box_zoom', 'reset', 'save'] , toolbar_location="above",
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

        #plot mem prediction
        if self.mempred != None:
            x_i = []
            x_m = []
            x_o = []

            for position, residue in enumerate(self.mempred):
                if residue.membrane_prediction == 'i':
                    #x_i.append(residue.x_i)
                    x_i.append(position + 1)
                elif residue.membrane_prediction == 'M':
                    x_m.append(position + 1)
                elif residue.membrane_prediction == 'o':
                    x_o.append(position + 1)

            y_i = x_i
            y_m = x_m
            y_o = x_o

            p2 = p.circle(x_i, y_i, size=10, color='green')
            p3 = p.circle(x_m, y_m, color='red', size=10)
            p4 = p.circle(x_o, y_o, color='yellow', size=10)

        #plot ss prediction
        if self.sspred != None:
            sequence = self.seq
            sequence_for_pred = []
            for indx, x in enumerate(sequence.seq):
                sequence_for_pred.append(x)

            xy_h = []
            xy_h_aa_name = []
            xy_c = []
            xy_c_aa_name = []
            xy_e = []
            xy_e_aa_name = []

            for position, residue in enumerate(self.sspred):
                if residue.ss2_prediction == 'H':
                    xy_h.append(position + 1)
                    for seq_position, aa in enumerate(sequence.seq):
                        if position == seq_position:
                            xy_h_aa_name.append(aa)
                elif residue.ss2_prediction == 'C':
                    xy_c.append(position + 1)
                    for seq_position, aa in enumerate(sequence.seq):
                        if position == seq_position:
                            xy_c_aa_name.append(aa)
                elif residue.ss2_prediction == 'E':  #todo CHECK!!
                    xy_e.append(position + 1)
                    for seq_position, aa in enumerate(sequence.seq):
                        if position == seq_position:
                            xy_e_aa_name.append(aa)

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
            y_dct = {}
            conserved_scores = []

            for x in self.conservationpred:
                conserved_scores.append(x.conservation_score)

            for score in consurf_values:
                temp = []
                for position, prediction in enumerate(conserved_scores):
                    if score == prediction:
                        temp.append(position+1)
                x_dct[score] = temp
                y_dct[score] = temp
                y_dct[score] = [int(x) + 3 for x in y_dct[score]]

            Blues = {1: '#f7fbff', 2: '#deebf7', 3: '#c6dbef', 4: '#9ecae1', 5: '#6baed6', 6: '#4292c6',
                         7: '#2171b5', 8: '#08519c', 9: '#08306b'}

            for x in Blues:
                p.triangle(x_dct[x], y_dct[x], size=5, color=Blues[x])

        # plot disorder prediction
        if self.dispred != None:

            x_d = []
            y_d = []
            x_nd = []
            y_nd = []

            for position, disorder_value in enumerate(self.dispred):
                if disorder_value.disorder_prediction >= 0.5:
                    x_d.append(position + 1)
                    y_d.append(position - 2)
                elif disorder_value.disorder_prediction < 0.5:
                    x_nd.append(position + 1)
                    y_nd.append(position - 2)

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



