from bokeh.plotting import figure, show, output_file
from bokeh.models import Range1d, Legend, ColumnDataSource
from bokeh.models.tools import HoverTool


class ContactPlotPlusFigure(object):

    def __init__(self, seq, conpred, factor=1):

        self.conpred = conpred
        self.conservpred = None
        self.mempred_i = None
        self.mempred_m = None
        self.mempred_o = None
        self.sspred_h = None
        self.sspred_c = None
        self.sspred_e = None
        self.dispred_o = None
        self.dispred_dis = None
        self.conservpred = None
        self.reference = None
        self.contact_w_reference_mis = None
        self.contact_w_reference_match = None
        self.reference_points = None
        self.contact_circles = None
        self.conpred.sequence = seq
        self.conpred.set_sequence_register()
        self.conpred.remove_neighbors(inplace=True)
        self.conpred.sort('raw_score', reverse=True, inplace=True)
        self.conpred = self.conpred[:self.conpred.sequence.seq_len * factor]
        self.canvas = figure(x_range=Range1d(0, self.conpred.ncontacts), y_range=Range1d(0, self.conpred.ncontacts),
                             plot_width=700, tools=['pan', 'wheel_zoom', 'box_zoom', 'reset', 'save'], plot_height=575,
                             toolbar_location="above")
        self.add_contact_layer()

    @property
    def spectrum(self):
        return {1: '#f7fbff',
                2: '#deebf7',
                3: '#c6dbef',
                4: '#9ecae1',
                5: '#6baed6',
                6: '#4292c6',
                7: '#2171b5',
                8: '#08519c',
                9: '#08306b'
                }

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
    def layers(self):
        return {"Mem_pred-I": [self.mempred_i],
                "Mem_pred-M": [self.mempred_m],
                "Mem_pred-O": [self.mempred_o],
                "SS-H": [self.sspred_h],
                "SS-C": [self.sspred_c],
                "SS-E": [self.sspred_e],
                "Disorder": [self.dispred_dis],
                "Ordered": [self.dispred_o],
                "Contact": [self.contact_circles],
                "Mismatch": [self.contact_w_reference_mis],
                "Match": [self.contact_w_reference_match],
                "Reference pdb": [self.reference_points]
                }

    @property
    def legend(self):

        # had to change conditional inequality otherwise legend unless all layers present
        return Legend(items=[(key, value) for key, value in self.layers.items() if value != [None]],
                      location='center')

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

    def add_reference_layer(self, reference):
        self.reference = reference
        id = [x.id for x in self.reference]
        res1 = [x.res1 for x in self.reference]
        res2 = [x.res2 for x in self.reference]
        res1_seq = [x.res1_seq for x in self.reference]
        res2_seq = [x.res2_seq for x in self.reference]
        source = ColumnDataSource(data=dict(id=id, res1=res1, res2=res2, res1_seq=res1_seq, res2_seq=res2_seq))

        # plot structural contacts
        self.reference_points = self.canvas.circle(x='res1_seq', y='res2_seq', size=3, color="#D8D6D6", source=source)
        mirror_image = self.canvas.circle(x='res2_seq', y='res1_seq', size=3, color="#D8D6D6", source=source)

        self.canvas.add_tools(
            HoverTool(renderers=[self.reference_points],
                      tooltips=[('Contact', '@id'), ('res1', '@res1'), ('res2', '@res2')]))

        self.canvas.add_tools(
            HoverTool(renderers=[mirror_image],
                      tooltips=[('Contact', '@id'), ('res1', '@res1'), ('res2', '@res2')]))

        # visualise matches/mismatches
        # could not get into a single line loop???
        if self.reference != None:
            res1_seq_match = []
            res2_seq_match = []
            res1_seq_mismatch = []
            res2_seq_mismatch = []
            for y in self.reference:
                for x in self.conpred:
                    if x.id == y.id:
                        res1_seq_match.append(x.res1_seq)
                        res2_seq_match.append(x.res2_seq)
                    elif x.id != y.id:
                        res1_seq_mismatch.append(x.res1_seq)
                        res2_seq_mismatch.append(x.res2_seq)

            self.contact_w_reference_mis = self.canvas.circle(res1_seq_mismatch, res2_seq_mismatch, size=3,
                                                              color='#DD4968')
            self.canvas.circle(res2_seq_mismatch, res1_seq_mismatch, size=3, color='#DD4968')

            self.contact_w_reference_match = self.canvas.circle(res1_seq_match, res2_seq_match, size=3, color='black')
            self.canvas.circle(res2_seq_match, res1_seq_match, size=3, color='black')

    def add_mempred_layer(self, mempred):

        self.mempred = mempred
        if len(mempred) == len(self.conpred.sequence.seq):
            self.mempred_i = self.canvas.circle([residue.id for residue in self.mempred.search('i')],
                                                [residue.id for residue in self.mempred.search('i')],
                                                size=10, color='green')
            self.mempred_m = self.canvas.circle([residue.id for residue in self.mempred.search('M')],
                                                [residue.id for residue in self.mempred.search('M')],
                                                color='red', size=10)
            self.mempred_o = self.canvas.circle([residue.id for residue in self.mempred.search('o')],
                                                [residue.id for residue in self.mempred.search('o')],
                                                color='yellow', size=10)
        else:
            raise ValueError('membrane prediction is not the same length as sequence')

    def add_sspred_layer(self, sspred):

        self.sspred = sspred
        if len(sspred) == len(self.conpred.sequence.seq):
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
        else:
            raise ValueError('ss prediction is not the same length as sequence')

    def add_dispred_layer(self, dispred):

        self.dispred = dispred

        if len(dispred) == len(self.conpred.sequence.seq):

            x_d = [residue.id for residue in self.dispred.search(0.5)]
            x_nd = [residue.id for residue in self.dispred if residue.id not in x_d]
            y_d = [residue.id for residue in self.dispred.search(0.5)]
            y_nd = [residue.id for residue in self.dispred if residue.id not in y_d]
            y_d = [int(y) - 3 for y in y_d]
            y_nd = [int(y) - 3 for y in y_nd]

            self.dispred_o = self.canvas.circle(x_nd, y_nd, color='grey', size=5)
            self.dispred_dis = self.canvas.circle(x_d, y_d, size=5, color='lime')
        else:
            raise ValueError('disorder prediction is not the same length as sequence')

    def add_conservpred_layer(self, conservpred):
        self.conservpred = conservpred

        if len(conservpred) == len(self.conpred.sequence.seq):
            consurf_values = list(range(1, 10))

            x_dct = {}
            y_dct = {}

            for score in consurf_values:
                x_dct[score] = [residue.id for residue in self.conservpred.search(score)]
                y_dct[score] = [int(x) + 3 for x in x_dct[score]]

            # for x in self.spectrum:
            #     self.canvas.triangle(x_dct[x], y_dct[x], size=5, color=self.spectrum[x])

            # couldn't get this into a single line loop????
            for x in self.spectrum.keys():
                self.canvas.triangle(x_dct[x], y_dct[x], size=5, color=self.spectrum[x])


        else:
            raise ValueError('conservation prediction is not the same length as sequence')

    def plot(self, fname):
        # print(self.layers)
        self.canvas.add_layout(self.legend, 'right')

        output_file(fname)
        show(self.canvas)
