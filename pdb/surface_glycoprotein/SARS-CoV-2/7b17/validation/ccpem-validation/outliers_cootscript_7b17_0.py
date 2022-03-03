
from __future__ import division
import cPickle
try :
  import gobject
except ImportError :
  gobject = None
import sys

dict_residue_prop_objects = {}
class coot_extension_gui (object) :
  def __init__ (self, title) :
    import gtk
    self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
    scrolled_win = gtk.ScrolledWindow()
    self.outside_vbox = gtk.VBox(False, 2)
    self.inside_vbox = gtk.VBox(False, 0)
    self.window.set_title(title)
    self.inside_vbox.set_border_width(0)
    self.window.add(self.outside_vbox)
    self.outside_vbox.pack_start(scrolled_win, True, True, 0)
    scrolled_win.add_with_viewport(self.inside_vbox)
    scrolled_win.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)

  def finish_window (self) :
    import gtk
    self.outside_vbox.set_border_width(2)
    ok_button = gtk.Button("  Close  ")
    self.outside_vbox.pack_end(ok_button, False, False, 0)
    ok_button.connect("clicked", lambda b: self.destroy_window())
    self.window.connect("delete_event", lambda a, b: self.destroy_window())
    self.window.show_all()

  def destroy_window (self, *args) :
    self.window.destroy()
    self.window = None

  def confirm_data (self, data) :
    for data_key in self.data_keys :
      outlier_list = data.get(data_key)
      if outlier_list is not None and len(outlier_list) > 0 :
        return True
    return False

  def create_property_lists (self, data) :
    import gtk
    for data_key in self.data_keys :
      outlier_list = data[data_key]
      if outlier_list is None or len(outlier_list) == 0 :
        continue
      else :
        frame = gtk.Frame(self.data_titles[data_key])
        vbox = gtk.VBox(False, 2)
        frame.set_border_width(6)
        frame.add(vbox)
        self.add_top_widgets(data_key, vbox)
        self.inside_vbox.pack_start(frame, False, False, 5)
        list_obj = residue_properties_list(
          columns=self.data_names[data_key],
          column_types=self.data_types[data_key],
          rows=outlier_list,
          box=vbox)
        ##save property list frame object
        dict_residue_prop_objects[data_key] = list_obj
# Molprobity result viewer
class coot_molprobity_todo_list_gui (coot_extension_gui) :
  data_keys = [ "clusters","rama", "rota", "cbeta", "probe", "smoc", "fdr",
               "fsc","diffmap","cablam",
               "jpred"]
  data_titles = { "clusters"  : "Outlier residue clusters",
                  "rama"  : "Ramachandran outliers",
                  "rota"  : "Rotamer outliers",
                  "cbeta" : "C-beta outliers",
                  "probe" : "Severe clashes",
                  "smoc"  : "Local density fit (SMOC)",
                  "fdr": "Backbone position score (FDR)",
                  "fsc": "Local density fit (FSC)",
                  "diffmap": "Model-map difference",
                  "cablam": "Ca geometry (CaBLAM)",
                  "jpred":"SS prediction"}
  data_names = { "clusters"  : ["Chain","Residue","Cluster","Outlier types"],
                 "rama"  : ["Chain", "Residue", "Name", "Score"],
                 "rota"  : ["Chain", "Residue", "Name", "Score"],
                 "cbeta" : ["Chain", "Residue", "Name", "Conf.", "Deviation"],
                 "probe" : ["Atom 1", "Atom 2", "Overlap"],
                 "smoc" : ["Chain", "Residue", "Name", "Score"],
                 "fdr" : ["Chain", "Residue", "Name", "Score"],
                 "fsc" : ["Chain", "Residue", "Name", "Score"],
                 "diffmap" : ["Chain", "Residue", "Name", "Score"],
                 "cablam" : ["Chain", "Residue","Name","recommendation","DSSP"],
                 "jpred" : ["Chain", "Residue","Name","predicted SS","current SS"]}
  if (gobject is not None) :
    data_types = {  "clusters" : [gobject.TYPE_STRING, gobject.TYPE_STRING,
                             gobject.TYPE_INT, gobject.TYPE_STRING,
                             gobject.TYPE_PYOBJECT],
                    "rama" : [gobject.TYPE_STRING, gobject.TYPE_STRING,
                             gobject.TYPE_STRING, gobject.TYPE_FLOAT,
                             gobject.TYPE_PYOBJECT],
                   "rota" : [gobject.TYPE_STRING, gobject.TYPE_STRING,
                             gobject.TYPE_STRING, gobject.TYPE_FLOAT,
                             gobject.TYPE_PYOBJECT],
                   "cbeta" : [gobject.TYPE_STRING, gobject.TYPE_STRING,
                              gobject.TYPE_STRING, gobject.TYPE_STRING,
                              gobject.TYPE_FLOAT, gobject.TYPE_PYOBJECT],
                   "probe" : [gobject.TYPE_STRING, gobject.TYPE_STRING,
                              gobject.TYPE_FLOAT, gobject.TYPE_PYOBJECT],
                   "smoc" : [gobject.TYPE_STRING, gobject.TYPE_STRING,
                              gobject.TYPE_STRING,gobject.TYPE_FLOAT,
                             gobject.TYPE_PYOBJECT],
                   "fdr" : [gobject.TYPE_STRING, gobject.TYPE_STRING,
                              gobject.TYPE_STRING,gobject.TYPE_FLOAT,
                             gobject.TYPE_PYOBJECT],
                   "fsc" : [gobject.TYPE_STRING, gobject.TYPE_STRING,
                              gobject.TYPE_STRING,gobject.TYPE_FLOAT,
                             gobject.TYPE_PYOBJECT],
                   "diffmap" : [gobject.TYPE_STRING, gobject.TYPE_STRING,
                              gobject.TYPE_STRING,gobject.TYPE_FLOAT,
                             gobject.TYPE_PYOBJECT],
                   "cablam" : [gobject.TYPE_STRING, gobject.TYPE_STRING,
                              gobject.TYPE_STRING,gobject.TYPE_STRING,
                             gobject.TYPE_STRING,gobject.TYPE_PYOBJECT],
                   "jpred" : [gobject.TYPE_STRING, gobject.TYPE_STRING,
                              gobject.TYPE_STRING,gobject.TYPE_STRING,
                             gobject.TYPE_STRING,gobject.TYPE_PYOBJECT]}
  else :
    data_types = dict([ (s, []) for s in ["clusters","rama","rota","cbeta","probe","smoc",
                                          "fdr","fsc","diffmap","cablam","jpred"] ])

  def __init__ (self, data_file=None, data=None) :
    assert ([data, data_file].count(None) == 1)
    if (data is None) :
      data = load_pkl(data_file)
    if not self.confirm_data(data) :
      return
    coot_extension_gui.__init__(self, "Validation To-do list")
    self.dots_btn = None
    self.dots2_btn = None
    self._overlaps_only = True
    self.window.set_default_size(420, 600)
    self.create_property_lists(data)
    self.finish_window()

  def add_top_widgets (self, data_key, box) :
    import gtk
    if data_key == "probe" :
      hbox = gtk.HBox(False, 2)
      self.dots_btn = gtk.CheckButton("Show Probe dots")
      hbox.pack_start(self.dots_btn, False, False, 5)
      self.dots_btn.connect("toggled", self.toggle_probe_dots)
      self.dots2_btn = gtk.CheckButton("Overlaps only")
      hbox.pack_start(self.dots2_btn, False, False, 5)
      self.dots2_btn.connect("toggled", self.toggle_all_probe_dots)
      self.dots2_btn.set_active(True)
      self.toggle_probe_dots()
      box.pack_start(hbox, False, False, 0)

  def toggle_probe_dots (self, *args) :
    if self.dots_btn is not None :
      show_dots = self.dots_btn.get_active()
      overlaps_only = self.dots2_btn.get_active()
      if show_dots :
        self.dots2_btn.set_sensitive(True)
      else :
        self.dots2_btn.set_sensitive(False)
      show_probe_dots(show_dots, overlaps_only)

  def toggle_all_probe_dots (self, *args) :
    if self.dots2_btn is not None :
      self._overlaps_only = self.dots2_btn.get_active()
      self.toggle_probe_dots()

class rsc_todo_list_gui (coot_extension_gui) :
  data_keys = ["by_res", "by_atom"]
  data_titles = ["Real-space correlation by residue",
                 "Real-space correlation by atom"]
  data_names = {}
  data_types = {}

class residue_properties_list (object) :
  def __init__ (self, columns, column_types, rows, box,
      default_size=(380,200)) :
    assert len(columns) == (len(column_types) - 1)
    if (len(rows) > 0) and (len(rows[0]) != len(column_types)) :
      raise RuntimeError("Wrong number of rows:\n%s" % str(rows[0]))
    import gtk
    ##adding a column type for checkbox (bool) before atom coordinate
    if gobject is not None:
        column_types = column_types[:-1]+[bool]+[column_types[-1]]
    
    self.liststore = gtk.ListStore(*column_types)
    self.listmodel = gtk.TreeModelSort(self.liststore)
    self.listctrl = gtk.TreeView(self.listmodel)
    self.listctrl.column = [None]*len(columns)
    self.listctrl.cell = [None]*len(columns)
    for i, column_label in enumerate(columns) :
      cell = gtk.CellRendererText()
      column = gtk.TreeViewColumn(column_label)
      self.listctrl.append_column(column)
      column.set_sort_column_id(i)
      column.pack_start(cell, True)
      column.set_attributes(cell, text=i)
    ##add a cell for checkbox
    cell1 = gtk.CellRendererToggle()
    cell1.connect ("toggled", self.on_selected_toggled)
    column = gtk.TreeViewColumn('Dealt with',cell1,active=i+1)
    self.listctrl.append_column(column)
    #column.set_sort_column_id(i+1)
    #column.pack_start(cell1, True)
    
    self.listctrl.get_selection().set_mode(gtk.SELECTION_SINGLE)
    for row in rows :
      row = row[:-1] + (False,)+(row[-1],)
      self.listmodel.get_model().append(row)
    self.listctrl.connect("cursor-changed", self.OnChange)
    sw = gtk.ScrolledWindow()
    w, h = default_size
    if len(rows) > 10 :
      sw.set_size_request(w, h)
    else :
      sw.set_size_request(w, 30 + (20 * len(rows)))
    sw.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
    box.pack_start(sw, False, False, 5)
    inside_vbox = gtk.VBox(False, 0)
    sw.add(self.listctrl)

  def OnChange (self, treeview) :
    import coot # import dependency
    selection = self.listctrl.get_selection()
    (model, tree_iter) = selection.get_selected()
    if tree_iter is not None :
      row = model[tree_iter]
      xyz = row[-1]
      if isinstance(xyz, tuple) and len(xyz) == 3 :
        set_rotation_centre(*xyz)
        set_zoom(30)
        graphics_draw()
  ##check box toggle
  def on_selected_toggled(self,renderer,path):
    if path is not None:
      model = self.listmodel.get_model()
      it = model.get_iter(path)
      #set toggle
      model[it][-2] = not model[it][-2]
      #set checkboxes for same residues in other lists
      try:
        chain = model[it][0]
        residue = model[it][1]
        for data_key in dict_residue_prop_objects:
          prop_obj = dict_residue_prop_objects[data_key]
          for row in prop_obj.listmodel.get_model():
            if data_key == 'probe':
              atom1_split = row[0].split()
              atom2_split = row[1].split()
              if atom1_split[0] == chain and atom1_split[1] == residue:
                row[-2] = model[it][-2]
              elif atom2_split[0] == chain and atom2_split[1] == residue:
                row[-2] = model[it][-2]
            elif row[0] == chain and row[1] == residue:
              row[-2] = model[it][-2]
      except IndexError: pass

  def check_chain_residue(self,chain,residue):
      pass
  
def show_probe_dots (show_dots, overlaps_only) :
  import coot # import dependency
  n_objects = number_of_generic_objects()
  sys.stdout.flush()
  if show_dots :
    for object_number in range(n_objects) :
      obj_name = generic_object_name(object_number)
      if overlaps_only and not obj_name in ["small overlap", "bad overlap"] :
        sys.stdout.flush()
        set_display_generic_object(object_number, 0)
      else :
        set_display_generic_object(object_number, 1)
  else :
    sys.stdout.flush()
    for object_number in range(n_objects) :
      set_display_generic_object(object_number, 0)

def load_pkl (file_name) :
  pkl = open(file_name, "rb")
  data = cPickle.load(pkl)
  pkl.close()
  return data
data = {}
data['rama'] = []
data['cbeta'] = []
data['fdr'] = []
data['fsc'] = []
data['diffmap'] = []
data['jpred'] = []
data['rota'] = [('A', ' 347 ', 'PHE', 0.17997401638022617, (313.446, 313.096, 318.563))]
data['clusters'] = [('A', '403', 1, 'side-chain clash', (306.78, 305.481, 306.801)), ('A', '404', 1, 'side-chain clash', (301.248, 305.191, 309.89)), ('A', '406', 1, 'side-chain clash', (306.78, 305.481, 306.801)), ('A', '407', 1, 'side-chain clash\nsmoc Outlier', (301.248, 305.191, 309.89)), ('A', '408', 1, 'smoc Outlier', (300.182, 300.699, 312.12)), ('A', '519', 2, 'cablam Outlier', (306.5, 291.2, 341.6)), ('A', '520', 2, 'smoc Outlier', (307.361, 294.297, 343.742)), ('A', '521', 2, 'smoc Outlier', (305.726, 297.456, 345.149)), ('A', '523', 2, 'smoc Outlier', (306.385, 303.225, 342.414)), ('A', '524', 2, 'side-chain clash', (303.093, 303.656, 338.346)), ('A', '453', 3, 'smoc Outlier', (315.685, 303.456, 309.546)), ('A', '491', 3, 'smoc Outlier', (322.112, 298.879, 305.923)), ('A', '492', 3, 'smoc Outlier', (320.445, 302.173, 306.899)), ('A', '493', 3, 'smoc Outlier', (318.167, 304.346, 304.769)), ('A', '347', 4, 'Rotamer', (313.446, 313.096, 318.563)), ('A', '401', 4, 'side-chain clash', (309.078, 312.35, 315.58)), ('A', '509', 4, 'side-chain clash', (309.078, 312.35, 315.58)), ('A', '338', 5, 'side-chain clash', (307.149, 313.296, 328.259)), ('A', '341', 5, 'side-chain clash', (307.149, 313.296, 328.259)), ('A', '601', 5, 'Planar group:C2:C7:C8:N2:O7', (302.205, 316.151, 326.173)), ('A', '335', 6, 'smoc Outlier', (304.691, 314.888, 338.448)), ('A', '361', 6, 'smoc Outlier', (306.461, 309.967, 340.335)), ('A', '362', 6, 'smoc Outlier', (302.765, 310.756, 339.94)), ('A', '454', 7, 'side-chain clash\nsmoc Outlier', (320.68, 295.148, 312.549)), ('A', '467', 7, 'side-chain clash', (320.68, 295.148, 312.549)), ('A', '498', 8, 'smoc Outlier', (309.953, 316.723, 302.619)), ('A', '499', 8, 'smoc Outlier', (306.983, 319.013, 303.381)), ('A', '462', 9, 'smoc Outlier', (313.53, 290.265, 319.154)), ('A', '464', 9, 'smoc Outlier', (312.483, 295.47, 322.521)), ('A', '364', 10, 'side-chain clash', (296.748, 314.72, 304.489)), ('A', '367', 10, 'side-chain clash', (296.748, 314.72, 304.489)), ('A', '413', 11, 'smoc Outlier', (301.692, 290.595, 315.721)), ('A', '414', 11, 'smoc Outlier', (302.699, 292.983, 312.924)), ('A', '444', 12, 'backbone clash', (315.28, 321.692, 303.054)), ('A', '445', 12, 'backbone clash', (315.28, 321.692, 303.054)), ('A', '429', 13, 'side-chain clash', (305.918, 298.602, 325.896)), ('A', '514', 13, 'side-chain clash', (305.918, 298.602, 325.896)), ('A', '374', 14, 'side-chain clash', (299.561, 312.802, 318.365)), ('A', '436', 14, 'side-chain clash', (299.561, 312.802, 318.365)), ('A', '345', 15, 'smoc Outlier', (314.017, 319.233, 322.151)), ('A', '459', 15, 'smoc Outlier', (316.876, 288.284, 310.743)), ('B', '22', 1, 'side-chain clash', (285.0, 311.072, 305.189)), ('B', '29', 1, 'side-chain clash', (282.214, 304.885, 307.612)), ('B', '34', 1, 'side-chain clash', (283.349, 305.523, 308.818)), ('B', '36', 1, 'side-chain clash', (285.0, 311.072, 305.189)), ('B', '73', 1, 'smoc Outlier', (277.209, 300.105, 306.909)), ('B', '78', 1, 'side-chain clash', (282.214, 304.885, 307.612)), ('B', '94', 1, 'side-chain clash', (287.693, 303.806, 307.495)), ('B', '187', 2, 'side-chain clash', (328.672, 305.87, 287.766)), ('B', '193', 2, 'side-chain clash', (328.672, 305.87, 287.766)), ('B', '204', 2, 'smoc Outlier', (328.746, 302.555, 281.231)), ('B', '205', 2, 'smoc Outlier', (326.405, 304.803, 283.233)), ('B', '206', 2, 'smoc Outlier', (326.973, 308.545, 283.554)), ('B', '217', 2, 'smoc Outlier', (325.905, 303.308, 277.047)), ('B', '100D', 3, 'side-chain clash\nbackbone clash', (294.727, 304.95, 313.553)), ('B', '100E', 3, 'backbone clash', (294.727, 304.95, 313.553)), ('B', '96', 3, 'side-chain clash', (294.488, 301.054, 312.654)), ('B', '229', 4, 'side-chain clash', (313.799, 303.13, 291.039)), ('B', '236', 4, 'smoc Outlier', (312.497, 301.912, 295.132)), ('B', '236L', 4, 'side-chain clash', (313.799, 303.13, 291.039)), ('B', '170', 5, 'side-chain clash', (319.809, 308.777, 286.452)), ('B', '214', 5, 'side-chain clash', (319.809, 308.777, 286.452)), ('B', '196', 6, 'side-chain clash', (325.389, 293.736, 285.63)), ('B', '199', 6, 'side-chain clash', (325.389, 293.736, 285.63)), ('B', '82', 7, 'side-chain clash', (281.631, 322.414, 309.536)), ('B', '82C', 7, 'side-chain clash', (281.631, 322.414, 309.536)), ('B', '60', 8, 'side-chain clash', (286.653, 318.562, 318.084)), ('B', '63', 8, 'side-chain clash', (286.653, 318.562, 318.084)), ('B', '111', 9, 'side-chain clash', (287.5, 328.792, 305.431)), ('B', '87', 9, 'side-chain clash', (287.5, 328.792, 305.431)), ('B', '19', 10, 'side-chain clash', (275.019, 315.741, 307.494)), ('B', '81', 10, 'side-chain clash', (275.019, 315.741, 307.494)), ('B', '140', 11, 'smoc Outlier', (310.925, 309.422, 284.777)), ('B', '160', 11, 'smoc Outlier', (313.74, 313.041, 285.432)), ('B', '188A', 12, 'side-chain clash', (325.892, 313.905, 291.208)), ('B', '207', 12, 'side-chain clash', (325.892, 313.905, 291.208)), ('B', '105', 13, 'cablam Outlier', (293.4, 312.2, 299.7)), ('B', '91', 13, 'side-chain clash', (296.748, 314.72, 304.489)), ('B', '145', 14, 'smoc Outlier', (315.233, 301.465, 273.207)), ('B', '146', 14, 'smoc Outlier', (314.729, 299.358, 270.119)), ('B', '5', 15, 'smoc Outlier', (287.499, 308.442, 298.766)), ('B', '241', 15, 'cablam Outlier', (307.6, 303.0, 282.2))]
data['probe'] = [(' B  34  MET  SD ', ' B  94  THR HG22', -0.665, (287.693, 303.806, 307.495)), (' B 188A SER  O  ', ' B 207  ARG  NH1', -0.66, (325.892, 313.905, 291.208)), (' A 454  ARG  NH1', ' A 467  ASP  OD2', -0.585, (320.68, 295.148, 312.549)), (' A 444  LYS  HD3', ' A 445  VAL  O  ', -0.577, (315.28, 321.692, 303.054)), (' B  19  ARG  NH1', ' B  81  GLN  OE1', -0.554, (275.019, 315.741, 307.494)), (' B  82  MET  HB3', ' B  82C LEU HD21', -0.551, (281.631, 322.414, 309.536)), (' B 196  SER  OG ', ' B 199  VAL HG22', -0.55, (325.389, 293.736, 285.63)), (' B  96  GLY  N  ', ' B 100D ASP  OD2', -0.521, (294.488, 301.054, 312.654)), (' A 403  ARG  HB3', ' A 406  GLU  OE1', -0.517, (306.78, 305.481, 306.801)), (' B  60  GLN  HB3', ' B  63  VAL HG22', -0.504, (286.653, 318.562, 318.084)), (' B 229  ALA  HB1', ' B 236L MET  HG2', -0.494, (313.799, 303.13, 291.039)), (' A 404  GLY  O  ', ' A 407  VAL HG22', -0.489, (301.248, 305.191, 309.89)), (' A 364  ASP  O  ', ' A 367  VAL HG12', -0.474, (298.314, 313.982, 332.195)), (' B 100D ASP  OD1', ' B 100E PHE  N  ', -0.468, (294.727, 304.95, 313.553)), (' A 503  VAL HG11', ' B  91  TYR  CE2', -0.467, (296.748, 314.72, 304.489)), (' B  34  MET  HG3', ' B  78  VAL HG11', -0.467, (283.349, 305.523, 308.818)), (' B 170  ILE HG21', ' B 214  VAL HG21', -0.465, (319.809, 308.777, 286.452)), (' A 524  VAL  O  ', ' A 524  VAL HG13', -0.449, (303.093, 303.656, 338.346)), (' A 429  PHE  HE1', ' A 514  SER  HG ', -0.449, (305.918, 298.602, 325.896)), (' B 187  ILE HD11', ' B 193  THR HG22', -0.448, (328.672, 305.87, 287.766)), (' A 338  PHE  O  ', ' A 341  VAL HG22', -0.444, (307.149, 313.296, 328.259)), (' A 374  PHE  HA ', ' A 436  TRP  HB3', -0.431, (299.561, 312.802, 318.365)), (' B  87  THR HG22', ' B 111  VAL  H  ', -0.428, (287.5, 328.792, 305.431)), (' B  22  CYS  HB3', ' B  78  VAL  CG2', -0.425, (284.043, 308.734, 305.669)), (' B  22  CYS  HB2', ' B  36  TRP  CZ2', -0.407, (285.0, 311.072, 305.189)), (' A 401  VAL HG22', ' A 509  ARG  HG2', -0.405, (309.078, 312.35, 315.58)), (' B  29  PHE  HZ ', ' B  78  VAL HG13', -0.4, (282.214, 304.885, 307.612))]
data['cablam'] = [('A', '519', 'HIS', 'check CA trace,carbonyls, peptide', 'bend\n--SS-', (306.5, 291.2, 341.6)), ('B', '105', 'GLN', 'check CA trace,carbonyls, peptide', ' \nE---E', (293.4, 312.2, 299.7)), ('B', '241', 'LYS', ' beta sheet', ' \n----E', (307.6, 303.0, 282.2))]
data['smoc'] = [('A', 335, u'LEU', -0.633, (304.691, 314.888, 338.448)), ('A', 345, u'THR', -0.39, (314.017, 319.233, 322.151)), ('A', 361, u'CYS', -0.582, (306.461, 309.967, 340.335)), ('A', 362, u'VAL', -0.527, (302.765, 310.756, 339.94)), ('A', 407, u'VAL', -0.341, (302.697, 303.559, 312.407)), ('A', 408, u'ARG', -0.37, (300.182, 300.699, 312.12)), ('A', 413, u'GLY', -0.465, (301.692, 290.595, 315.721)), ('A', 414, u'GLN', -0.5, (302.699, 292.983, 312.924)), ('A', 453, u'TYR', -0.138, (315.685, 303.456, 309.546)), ('A', 454, u'ARG', -0.2, (317.123, 299.979, 309.06)), ('A', 459, u'SER', -0.092, (316.876, 288.284, 310.743)), ('A', 462, u'LYS', -0.051, (313.53, 290.265, 319.154)), ('A', 464, u'PHE', -0.218, (312.483, 295.47, 322.521)), ('A', 491, u'PRO', -0.362, (322.112, 298.879, 305.923)), ('A', 492, u'LEU', -0.356, (320.445, 302.173, 306.899)), ('A', 493, u'GLN', -0.328, (318.167, 304.346, 304.769)), ('A', 498, u'GLN', -0.106, (309.953, 316.723, 302.619)), ('A', 499, u'PRO', -0.296, (306.983, 319.013, 303.381)), ('A', 520, u'ALA', -0.507, (307.361, 294.297, 343.742)), ('A', 521, u'PRO', -1.0, (305.726, 297.456, 345.149)), ('A', 523, u'THR', -0.637, (306.385, 303.225, 342.414)), ('A', 335, u'LEU', -0.633, (304.691, 314.888, 338.448)), ('A', 345, u'THR', -0.39, (314.017, 319.233, 322.151)), ('A', 361, u'CYS', -0.582, (306.461, 309.967, 340.335)), ('A', 362, u'VAL', -0.527, (302.765, 310.756, 339.94)), ('A', 407, u'VAL', -0.341, (302.697, 303.559, 312.407)), ('A', 408, u'ARG', -0.37, (300.182, 300.699, 312.12)), ('A', 413, u'GLY', -0.465, (301.692, 290.595, 315.721)), ('A', 414, u'GLN', -0.5, (302.699, 292.983, 312.924)), ('A', 453, u'TYR', -0.138, (315.685, 303.456, 309.546)), ('A', 454, u'ARG', -0.2, (317.123, 299.979, 309.06)), ('A', 459, u'SER', -0.092, (316.876, 288.284, 310.743)), ('A', 462, u'LYS', -0.051, (313.53, 290.265, 319.154)), ('A', 464, u'PHE', -0.218, (312.483, 295.47, 322.521)), ('A', 491, u'PRO', -0.362, (322.112, 298.879, 305.923)), ('A', 492, u'LEU', -0.356, (320.445, 302.173, 306.899)), ('A', 493, u'GLN', -0.328, (318.167, 304.346, 304.769)), ('A', 498, u'GLN', -0.106, (309.953, 316.723, 302.619)), ('A', 499, u'PRO', -0.296, (306.983, 319.013, 303.381)), ('A', 520, u'ALA', -0.507, (307.361, 294.297, 343.742)), ('A', 521, u'PRO', -1.0, (305.726, 297.456, 345.149)), ('A', 523, u'THR', -0.637, (306.385, 303.225, 342.414)), ('A', 335, u'LEU', -0.633, (304.691, 314.888, 338.448)), ('A', 345, u'THR', -0.39, (314.017, 319.233, 322.151)), ('A', 361, u'CYS', -0.582, (306.461, 309.967, 340.335)), ('A', 362, u'VAL', -0.527, (302.765, 310.756, 339.94)), ('A', 407, u'VAL', -0.341, (302.697, 303.559, 312.407)), ('A', 408, u'ARG', -0.37, (300.182, 300.699, 312.12)), ('A', 413, u'GLY', -0.465, (301.692, 290.595, 315.721)), ('A', 414, u'GLN', -0.5, (302.699, 292.983, 312.924)), ('A', 453, u'TYR', -0.138, (315.685, 303.456, 309.546)), ('A', 454, u'ARG', -0.2, (317.123, 299.979, 309.06)), ('A', 459, u'SER', -0.092, (316.876, 288.284, 310.743)), ('A', 462, u'LYS', -0.051, (313.53, 290.265, 319.154)), ('A', 464, u'PHE', -0.218, (312.483, 295.47, 322.521)), ('A', 491, u'PRO', -0.362, (322.112, 298.879, 305.923)), ('A', 492, u'LEU', -0.356, (320.445, 302.173, 306.899)), ('A', 493, u'GLN', -0.328, (318.167, 304.346, 304.769)), ('A', 498, u'GLN', -0.106, (309.953, 316.723, 302.619)), ('A', 499, u'PRO', -0.296, (306.983, 319.013, 303.381)), ('A', 520, u'ALA', -0.507, (307.361, 294.297, 343.742)), ('A', 521, u'PRO', -1.0, (305.726, 297.456, 345.149)), ('A', 523, u'THR', -0.637, (306.385, 303.225, 342.414)), ('B', 5, u'VAL', -0.688, (287.499, 308.442, 298.766)), ('B', 73, u'ASN', -0.507, (277.209, 300.105, 306.909)), ('B', 140, u'LEU', -0.296, (310.925, 309.422, 284.777)), ('B', 145, u'GLY', -0.381, (315.233, 301.465, 273.207)), ('B', 146, u'GLY', -0.473, (314.729, 299.358, 270.119)), ('B', 160, u'ALA', -0.212, (313.74, 313.041, 285.432)), ('B', 204, u'THR', -0.547, (328.746, 302.555, 281.231)), ('B', 205, u'ILE', -0.415, (326.405, 304.803, 283.233)), ('B', 206, u'SER', -0.474, (326.973, 308.545, 283.554)), ('B', 217, u'GLN', -0.473, (325.905, 303.308, 277.047)), ('B', 236, u'MET', -0.063, (312.497, 301.912, 295.132)), ('B', 5, u'VAL', -0.688, (287.499, 308.442, 298.766)), ('B', 73, u'ASN', -0.507, (277.209, 300.105, 306.909)), ('B', 140, u'LEU', -0.296, (310.925, 309.422, 284.777)), ('B', 145, u'GLY', -0.381, (315.233, 301.465, 273.207)), ('B', 146, u'GLY', -0.473, (314.729, 299.358, 270.119)), ('B', 160, u'ALA', -0.212, (313.74, 313.041, 285.432)), ('B', 204, u'THR', -0.547, (328.746, 302.555, 281.231)), ('B', 205, u'ILE', -0.415, (326.405, 304.803, 283.233)), ('B', 206, u'SER', -0.474, (326.973, 308.545, 283.554)), ('B', 217, u'GLN', -0.473, (325.905, 303.308, 277.047)), ('B', 236, u'MET', -0.063, (312.497, 301.912, 295.132)), ('B', 5, u'VAL', -0.688, (287.499, 308.442, 298.766)), ('B', 73, u'ASN', -0.507, (277.209, 300.105, 306.909)), ('B', 140, u'LEU', -0.296, (310.925, 309.422, 284.777)), ('B', 145, u'GLY', -0.381, (315.233, 301.465, 273.207)), ('B', 146, u'GLY', -0.473, (314.729, 299.358, 270.119)), ('B', 160, u'ALA', -0.212, (313.74, 313.041, 285.432)), ('B', 204, u'THR', -0.547, (328.746, 302.555, 281.231)), ('B', 205, u'ILE', -0.415, (326.405, 304.803, 283.233)), ('B', 206, u'SER', -0.474, (326.973, 308.545, 283.554)), ('B', 217, u'GLN', -0.473, (325.905, 303.308, 277.047)), ('B', 236, u'MET', -0.063, (312.497, 301.912, 295.132))]
handle_read_draw_probe_dots_unformatted("/home/ccpem/agnel/gisaid/countries_seq/structure_data/emdb/EMD-11980/7b17/Model_validation_1/validation_cootdata/molprobity_probe7b17_0.txt", 0, 0)
show_probe_dots(True, True)
gui = coot_molprobity_todo_list_gui(data=data)
