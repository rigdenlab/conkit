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

    @property
    def aa_properites(self):
        return {'A': 'Non-polar, aliphatic',
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

        #plot mem prediction
        if self.mempred != None:
            p2 = p.circle([residue.id for residue in self.mempred.search('i')],
                          [residue.id for residue in self.mempred.search('i')], size=10, color='green')
            p3 = p.circle([residue.id for residue in self.mempred.search('M')],
                          [residue.id for residue in self.mempred.search('M')], color='red', size=10)
            p4 = p.circle([residue.id for residue in self.mempred.search('o')],
                          [residue.id for residue in self.mempred.search('o')], color='yellow', size=10)


        #plot ss prediction
        if self.sspred != None:
            sequence = self.seq
            xy_h = [residue.id for residue in self.sspred.search('H')]
            xy_c = [residue.id for residue in self.sspred.search('C')]
            xy_e = [residue.id for residue in self.sspred.search('E')]
            xy_h_aa_name = []
            xy_c_aa_name = []
            xy_e_aa_name = []

            for id in xy_h:
                for seq_position, aa in enumerate(sequence.seq):
                    if int(id) == int(seq_position)+1:
                        xy_h_aa_name.append(aa)
            for id in xy_c:
                for seq_position, aa in enumerate(sequence.seq):
                    if int(id) == int(seq_position)+1:
                        xy_c_aa_name.append(aa)
            for id in xy_e:
                for seq_position, aa in enumerate(sequence.seq):
                    if int(id) == int(seq_position)+1:
                        xy_e_aa_name.append(aa)

            #add aa properties to each residue from aa property dictionary
            h_aa_properies = []

            for x in xy_h_aa_name:
                for y in self.aa_properites:
                    if x == y:
                        h_aa_properies.append(self.aa_properites[y])

            c_aa_properies = []
            for x in xy_c_aa_name:
                for y in self.aa_properites:
                    if x == y:
                        c_aa_properies.append(self.aa_properites[y])

            e_aa_properies = []
            for x in xy_e_aa_name:
                for y in self.aa_properites:
                    if x == y:
                        e_aa_properies.append(self.aa_properites[y])

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
            conserved_scores = [residue.prediction for residue in self.conservationpred.search()]

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

            disorder = [residue.prediction for residue in self.dispred.search()]

            for position, prediction in enumerate(disorder):
                if prediction >= 0.5:
                    x_d.append(position + 1)
                    y_d.append(position - 2)
                elif prediction < 0.5:
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

if __name__ == '__main__':

    import conkit.io

    conpred = conkit.io.read('/Users/shahrammesdaghi/Downloads/w9dy28_meta_respre.evfold', 'evfold').top
    seq = conkit.io.read('/Users/shahrammesdaghi/Downloads/w9dy28.fasta', 'fasta').top
    sspred = conkit.io.read('/Users/shahrammesdaghi/Downloads/w9dy28.ss2', 'psipred')
    conservpred = conkit.io.read('/Users/shahrammesdaghi/Downloads/consurf.grades.txt', 'consurf')
    mempred = conkit.io.read('/Users/shahrammesdaghi/Downloads/query.result.txt', 'topcons')
    dispred = conkit.io.read('/Users/shahrammesdaghi/Downloads/iupred2a.txt', 'iupred2a')

    #test_plot = ContactPlotPlusFigure(seq, conpred, ss_prediction, con_pred, mem_pred, dis_pred)
    #test_plot.plot_plus()

    test_plot_2 = conkit.plot.ContactPlotPlusFigure(seq=seq, conpred=conpred, conservationpred=conservpred,
                                                    mempred=mempred, dispred=dispred, sspred=sspred)


    test_plot_2.plot()

