from conkit.plot.figure import Figure
from bokeh.plotting import figure, show, output_file
from bokeh.models import Range1d, Legend, ColumnDataSource
from bokeh.models.tools import HoverTool


class ContactPlotPlusFigure(Figure):

    def __init__(self, seq, conpred, factor=1):

        self._conpred = conpred
        self._conservpred = None
        self._mempred_i = None
        self._mempred_m = None
        self._mempred_o = None
        self._sspred_h = None
        self._sspred_c = None
        self._sspred_e = None
        self._dispred_o = None
        self._dispred_dis = None
        self._conservpred = None
        self._contact_circles = None
        self.conpred.sequence = seq
        self.conpred.set_sequence_register()
        self.conpred.remove_neighbors(inplace=True)
        self.conpred.sort('raw_score', reverse=True, inplace=True)
        self.conpred = self.conpred[:self.conpred.sequence.seq_len * factor]
        self._canvas = figure(x_range=Range1d(0, self.conpred.ncontacts), y_range=Range1d(0, self.conpred.ncontacts),
                              plot_width=700, tools=['pan', 'wheel_zoom', 'box_zoom', 'reset', 'save'], plot_height=575,
                              toolbar_location="above")
        self.add_contact_layer()

    @property
    def conpred(self):
        return self._conpred

    @conpred.setter
    def conpred(self, value):
        self._conpred = value

    @property
    def canvas(self):
        return self._canvas

    @canvas.setter
    def canvas(self, value):
        self._canvas = value

    @property
    def dispred_o(self):
        return self._dispred_o

    @dispred_o.setter
    def dispred_o(self, value):
        self._dispred_o = value

    @property
    def dispred_dis(self):
        return self._dispred_dis

    @dispred_dis.setter
    def dispred_dis(self, value):
        self._dispred_dis = value

    @property
    def mempred_i(self):
        return self._mempred_i

    @mempred_i.setter
    def mempred_i(self, value):
        self._mempred_i = value

    @property
    def mempred_m(self):
        return self._mempred_m

    @mempred_m.setter
    def mempred_m(self, value):
        self._mempred_m = value

    @property
    def mempred_o(self):
        return self._mempred_o

    @mempred_o.setter
    def mempred_o(self, value):
        self._mempred_o = value

    @property
    def sspred_h(self):
        return self._sspred_h

    @sspred_h.setter
    def sspred_h(self, value):
        self._sspred_h = value

    @property
    def sspred_c(self):
        return self._sspred_c

    @sspred_c.setter
    def sspred_c(self, value):
        self._sspred_c = value

    @property
    def sspred_e(self):
        return self._sspred_e

    @sspred_e.setter
    def sspred_e(self, value):
        self._sspred_e = value

    @property
    def conservpred(self):
        return self._conservpred

    @conservpred.setter
    def conservpred(self, value):
        self._conservpred = value

    @property
    def contact_circles(self):
        return self._contact_circles

    @contact_circles.setter
    def contact_circles(self, value):
        self._contact_circles = value

    @property
    def spectrum(self):
        return {1: '#f7fbff', 2: '#deebf7', 3: '#c6dbef', 4: '#9ecae1', 5: '#6baed6', 6: '#4292c6', 7: '#2171b5',
                8: '#08519c', 9: '#08306b'}

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

    @property
    def legend(self):

        # construct legend
        legend_data = {"Mem_pred-I": [self.mempred_i], "Mem_pred-M": [self.mempred_m], "Mem_pred-O": [self.mempred_o],
                       "SS-H": [self.sspred_h], "SS-C": [self.sspred_c], "SS-E": [self.sspred_e],
                       "Disorder": [self.dispred_dis],
                       "Ordered": [self.dispred_o]}

        items = [("Contact", [self.contact_circles])]

        for key, value in legend_data.items():
            if value != [None]:
                legend_item = (key, value)
                items.append(legend_item)

        return Legend(items=items, location="center")

    def add_contact_layer(self):

        id = [x.id for x in self.conpred]
        res1 = [x.res1 for x in self.conpred]
        res1_chain = [x.res1_chain for x in self.conpred]
        res1_seq = [x.res1_seq for x in self.conpred]
        res2 = [x.res2 for x in self.conpred]
        res2_chain = [x.res2_chain for x in self.conpred]
        res2_seq = [x.res2_seq for x in self.conpred]
        raw_score = [x.raw_score for x in self.conpred]

        source = ColumnDataSource(data=dict(id=id, res1=res1, res1_chain=res1_chain, res1_seq=res1_seq, res2=res2,
                                            res2_chain=res2_chain, res2_seq=res2_seq, raw_score=raw_score))

        self.contact_circles = self.canvas.circle(x='res1_seq', y='res2_seq', size=3, color='black', source=source)

        self.canvas.add_tools(HoverTool(renderers=[self.contact_circles],
                                        tooltips=[('Contact', '@id'), ('res1', '@res1'), ('res2', '@res2'),
                                                  ('confidence', '@raw_score')]))

        p1a = self.canvas.circle(x='res2_seq', y='res1_seq', size=3, color='black', source=source)

        self.canvas.add_tools(
            HoverTool(renderers=[p1a], tooltips=[('Contact', '@id'), ('res1', '@res1'), ('res2', '@res2'),
                                                 ('confidence', '@raw_score')]))

    def add_mempred_layer(self, mempred):

        self.mempred = mempred
        self.mempred_i = self.canvas.circle([residue.id for residue in self.mempred.search('i')],
                                            [residue.id for residue in self.mempred.search('i')], size=10,
                                            color='green')
        self.mempred_m = self.canvas.circle([residue.id for residue in self.mempred.search('M')],
                                            [residue.id for residue in self.mempred.search('M')], color='red', size=10)
        self.mempred_o = self.canvas.circle([residue.id for residue in self.mempred.search('o')],
                                            [residue.id for residue in self.mempred.search('o')], color='yellow',
                                            size=10)

    def add_sspred_layer(self, sspred):

        self.sspred = sspred
        # find residue positions of specific ss prediction
        xy_h = [int(residue.id) for residue in self.sspred.search('H')]
        xy_c = [int(residue.id) for residue in self.sspred.search('C')]
        xy_e = [int(residue.id) for residue in self.sspred.search('E')]

        # find the residue name for those specific positions
        xy_h_aa_name = [self.conpred.sequence.seq[res_idx - 1] for res_idx in xy_h]
        xy_c_aa_name = [self.conpred.sequence.seq[res_idx - 1] for res_idx in xy_c]
        xy_e_aa_name = [self.conpred.sequence.seq[res_idx - 1] for res_idx in xy_e]

        # add aa properties to each residue from aa property dictionary
        h_aa_properies = [self.aa_properites[x] for x in xy_h_aa_name]
        c_aa_properies = [self.aa_properites[x] for x in xy_c_aa_name]
        e_aa_properies = [self.aa_properites[x] for x in xy_e_aa_name]

        # plot ss on to the contact plot with residue properties
        source2 = ColumnDataSource(data=dict(x=xy_h, y=xy_h, aa_name=xy_h_aa_name, aa_properties=h_aa_properies))
        source3 = ColumnDataSource(data=dict(x=xy_c, y=xy_c, aa_name=xy_c_aa_name, aa_properties=c_aa_properies))
        source4 = ColumnDataSource(data=dict(x=xy_e, y=xy_e, aa_name=xy_e_aa_name, aa_properties=e_aa_properies))

        self.sspred_h = self.canvas.circle(x='x', y='y', size=3, color='orange', source=source2)
        self.canvas.add_tools(
            HoverTool(renderers=[self.sspred_h], tooltips=[('Position', '@x'), ('Residue Name', '@aa_name'),
                                                           ('Properties', '@aa_properties')]))
        self.sspred_c = self.canvas.circle(x='x', y='y', color='blue', size=3, source=source3)
        self.canvas.add_tools(
            HoverTool(renderers=[self.sspred_c], tooltips=[('Position', '@x'), ('Residue Name', '@aa_name'),
                                                           ('Properties', '@aa_properties')]))
        self.sspred_e = self.canvas.circle(x='x', y='y', color='pink', size=3, source=source4)
        self.canvas.add_tools(
            HoverTool(renderers=[self.sspred_e], tooltips=[('Position', '@x'), ('Residue Name', '@aa_name'),
                                                           ('Properties', '@aa_properties')]))

    def add_dispred_layer(self, dispred):

        self.dispred = dispred

        # by default plot of 'ordered'
        x_nd = [i for i in range(1, self.conpred.sequence.seq_len + 1)]
        y_nd = [i for i in range(-2, self.conpred.sequence.seq_len - 2)]
        self.dispred_o = self.canvas.circle(x_nd, y_nd, color='grey', size=5)

        # overwirites 'ordered' residues with disordered residues
        x_d = [residue.id for residue in self.dispred.search(0.5)]
        y_d = [(int(residue.id) - 3) for residue in self.dispred.search(0.5)]
        self.dispred_dis = self.canvas.circle(x_d, y_d, size=5, color='lime')

    def add_conservpred_layer(self, conservpred):
        self.conservpred = conservpred
        consurf_values = list(range(1, 10))

        x_dct = {}
        y_dct = {}

        for score in consurf_values:
            x_dct[score] = [residue.id for residue in self.conservpred.search(score)]
            y_dct[score] = [int(x) + 3 for x in x_dct[score]]

        for x in self.spectrum:
            self.canvas.triangle(x_dct[x], y_dct[x], size=5, color=self.spectrum[x])

    def plot(self, fname):
        self.add_contact_layer()
        self.canvas.add_layout(self.legend, 'right')
        output_file(fname)
        show(self.canvas)


if __name__ == '__main__':
    import conkit.io

    conpred = conkit.io.read('/Users/shahrammesdaghi/Downloads/w9dy28_meta_respre.evfold', 'evfold').top
    seq = conkit.io.read('/Users/shahrammesdaghi/Downloads/w9dy28.fasta', 'fasta').top
    sspred = conkit.io.read('/Users/shahrammesdaghi/Downloads/w9dy28.ss2', 'psipred')
    conservpred = conkit.io.read('/Users/shahrammesdaghi/Downloads/consurf.grades.txt', 'consurf')
    mempred = conkit.io.read('/Users/shahrammesdaghi/Downloads/query.result.txt', 'topcons')
    dispred = conkit.io.read('/Users/shahrammesdaghi/Downloads/iupred2a.txt', 'iupred2a')

    test_plot_2 = conkit.plot.ContactPlotPlusFigure(seq=seq, conpred=conpred)

    test_plot_2.add_conservpred_layer(conservpred)
    test_plot_2.add_dispred_layer(dispred)
    test_plot_2.add_mempred_layer(mempred)
    test_plot_2.add_sspred_layer(sspred)
    test_plot_2.plot(fname='/Users/shahrammesdaghi/Downloads/test.html')
