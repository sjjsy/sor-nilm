#!/usr/bin/env python3
"""analyze.py -- An energy consumption data analyzer

Usage: analyze.py ACTION ARGS

Actions:
  convert CSV_FILES   Converts the given CSV files to DB files by device.
  plot DB_FILES       Plots a graph for each DB file depicting its content.
  hmm                 Execute our HMM test.

Examples:
  $ src/analyze.py convert in/mdr_2016-03_export_9*.csv
  $ src/analyze.py plot db/D22.db
  $ src/analyze.py hmm
"""

import sys

from pathlib import Path
from datetime import datetime, timedelta

def perr(msg):
  print('ERR: %s' % (msg))

def abort():
  perr('Critical failure! Aborting...')
  sys.exit(2)

class Row:
  def __init__(self, pf, irow, device_id, type_id, value_time_stamp, data_value):
    self.pf = pf
    self.irow = irow
    self.device_id = device_id
    self.type_id = type_id
    self.value_time_stamp = value_time_stamp
    self.data_value = data_value
    self.dt = datetime.strptime(value_time_stamp, '%Y-%m-%d %H:%M:%S')

  def __str__(self):
    return '%s, %s, %s, %s, %s, %s' % (self.pf, self.irow, self.device_id, self.type_id, self.value_time_stamp, self.data_value)

class Measurement:
  def __init__(self, dt, p1=None, p2=None, p3=None):
    self.dt = dt
    self.phases = [p1, p2, p3]
  @property
  def p1(self):
    return self.phases[0]
  @property
  def p2(self):
    return self.phases[1]
  @property
  def p3(self):
    return self.phases[2]
  @property
  def pwr(self):
    return sum(self.phases) * 230
  def __str__(self):
    return '[%s][%s+%s+%s]' % (self.dt.strftime('%Y-%m-%d %H:%M:%S'), self.p1, self.p2, self.p3)
  def to_csv(self):
    return '%s;%.1f;%.1f;%.1f' % (self.dt.strftime('%Y-%m-%d %H:%M:%S'), self.p1, self.p2, self.p3)
  @classmethod
  def from_csv(self, csvrow):
    d, a, b, c = csvrow.split(';')
    return Measurement(datetime.strptime(d, '%Y-%m-%d %H:%M:%S'), float(a), float(b), float(c))

class MeasurementSequence:
  def __init__(self, name, pfdb=None, loaddb=True):
    self.name = name
    self.vmmnt = []
    self.dmmnt = {}
    self.pfdb = pfdb
    if loaddb and self.pfdb and self.pfdb.exists():
      self.load()
  def __str__(self):
    return '%s: [%s - %s], N=%d' % (self.name, self.vmmnt[0].dt.strftime('%m-%d %H:%M:%S'), self.vmmnt[-1].dt.strftime('%m-%d %H:%M:%S'), len(self.vmmnt))
  @property
  def n(self):
    return len(self.vmmnt)
  @property
  def td(self):
    return self.vmmnt[-1].dt - self.vmmnt[0].dt if self.n > 1 else timedelta(seconds=0)
  def addm(self, mmnt):
    if mmnt.dt not in self.dmmnt:
      self.vmmnt.append(mmnt)
      self.dmmnt[mmnt.dt] = mmnt
  def save(self):
    print('Saving DB %s to %s...' % (self, self.pfdb))
    with self.pfdb.open('w') as f:
      for mmnt in self.vmmnt:
        f.write(mmnt.to_csv() + '\n')
  def load(self):
    print('Loading DB at %s...' % (self.pfdb))
    with self.pfdb.open('r') as f:
      for ln in f.readlines():
        mmnt = Measurement.from_csv(ln)
        self.vmmnt.append(mmnt)
        self.dmmnt[mmnt.dt] = mmnt
    print('Loaded DB %s!' % (self))
  def compstats(self):
    if self.n == 0:
      self.dt_avg = None
      self.p1_avg = self.p2_avg = self.p3_avg = 0
      self.mmnt_avg = Measurement(self.dt_avg, self.p1_avg, self.p2_avg, self.p3_avg)
      return
    mmnt1 = self.vmmnt[0]
    td = timedelta(seconds=0)
    p1 = p2 = p3 = 0
    for mmnt in self.vmmnt:
      td += mmnt.dt - mmnt1.dt
      p1 += mmnt.p1
      p2 += mmnt.p2
      p3 += mmnt.p3
    n = len(self.vmmnt)
    self.dt_avg = mmnt1.dt + td/n
    self.p1_avg = p1/n
    self.p2_avg = p2/n
    self.p3_avg = p3/n
    self.mmnt_avg = Measurement(self.dt_avg, self.p1_avg, self.p2_avg, self.p3_avg)
  def subseqtd(self, td=60):
    pfdb = self.pfdb.parent / (self.pfdb.stem + '-TD%d.dbv' % (td))
    mstd = MeasurementSequence(pfdb.stem, pfdb, loaddb=False)
    mtmp = self.vmmnt[0]
    nmstmp = 'temp-%d' % (td)
    mstmp = MeasurementSequence(nmstmp)
    td = timedelta(seconds=td)
    #print(td)
    for mmnt in self.vmmnt:
      mstmp.addm(mmnt)
      if mmnt.dt - mtmp.dt >= td:
        mtmp = mmnt
        mstmp.compstats()
        #print(mstmp, mstmp.mmnt_avg)
        mstd.addm(mstmp.mmnt_avg)
        mstmp = MeasurementSequence(nmstmp)
    if mstmp.n > 0:
      mstmp.compstats()
      mstd.addm(mstmp.mmnt_avg)
    print(mstd)
    return mstd

class Device(MeasurementSequence):
  def __init__(self, iid):
    self.iid = iid
    super().__init__('D%d' % (self.iid), Path('db/D%s.db' % (self.iid)))

action = sys.argv[1] if len(sys.argv) > 0 else None

### PARSING
if action == 'convert':
  ddevice = {}
  for arg in sys.argv[2:]:
    ## Parse the datafiles
    pi = Path(arg)
    if not pi.exists():
      perr('Cannot convert nonexisting file %s!' % (pi))
      continue
    #rows = []
    print('Parsing file %s...' % (pi))
    irow = irowok = 0
    with pi.open('r') as f:
      s = f.read()
      for row in s.split('\n'):
        irow += 1
        # Parsing the columns of the data row and converting it into a Row object
        cols = row.split(';')
        if len(cols) < 4:
          if len(row) > 2:
            perr('[R%d]: Invalid row: %s!' % (irow, row))
          continue
        if irow == 1 and cols[0] not in '1234567890':
          print('[R%d]: Skipping header %s...' % (irow, row))
          continue
        device_id, type_id, value_time_stamp, data_value = cols
        row = Row(pi, irow, int(device_id), int(type_id), value_time_stamp, float(data_value))
        # Find the device which this measurement represents
        if row.device_id not in ddevice:
          device = ddevice[row.device_id] = Device(row.device_id)
        else:
          device = ddevice[row.device_id]
        # Find the measurement which this row represents
        if row.dt not in device.dmmnt:
          measurement = Measurement(row.dt)
        else:
          measurement = device.dmmnt[row.dt]
        # Complement the measurement with the phase value
        phase = row.type_id - 1
        if measurement.phases[phase] != None:
          #perr('[R%d]: Found a duplicate phase data value: %s for P%d at %s!' % (irow, row.data_value, row.type_id, measurement))
          continue
        measurement.phases[phase] = row.data_value
        # Store the mesurement
        device.addm(measurement)
        irowok += 1
      print('Processed %d (-%d) rows from %s...' % (irow, irow - irowok, pi))
  ## Prune the device objects
  for device in ddevice.values():
    vmmnt = []
    for mmnt in device.vmmnt:
      if None in mmnt.phases:
        #perr('Missing phase(s) from %s...' % (mmnt))
        continue
      vmmnt.append(mmnt)
    device.vmmnt = sorted(vmmnt, key=lambda m: m.dt)
    device.dmmnt = {mmnt.dt: mmnt for mmnt in device.vmmnt}
  ## Save the device sequences and produce and save their avg variants
  for device in ddevice.values():
    device.save()
    mssq = device.subseqtd(60)
    mssq.save()
    mssq = device.subseqtd(3600)
    mssq.save()

### DB VARIANTS
elif action == 'dbvar':
  for arg in sys.argv[2:]:
    pfdb = Path(arg)
    if not pfdb.exists():
      perr('No DB at %s!' % (pfdb))
      continue
    mseq = MeasurementSequence(pfdb.stem, pfdb)
    mssq = mseq.subseqtd(60)
    mssq.save()
    mssq = mseq.subseqtd(3600)
    mssq.save()

### PLOTTING
elif action == 'plot':
  import matplotlib.pyplot as plt
  import matplotlib.dates as mdates
  import matplotlib.ticker as mtckr
  for arg in sys.argv[2:]:
    pfdb = Path(arg)
    if not pfdb.exists():
      perr('No DB at %s!' % (pfdb))
      continue
    mseq = MeasurementSequence(pfdb.stem, pfdb)
    print('Plotting %s...' % (mseq))
    fig, ax = plt.subplots(figsize=(16,9))
    mseq.compstats()
    ax.set_title(str(mseq) + ', AVG:%d+%d+%d=%d' % (mseq.p1_avg, mseq.p2_avg, mseq.p3_avg, mseq.p1_avg+mseq.p2_avg+mseq.p3_avg))
    ax.plot([mmnt.dt for mmnt in mseq.vmmnt], [mmnt.p1 for mmnt in mseq.vmmnt], label='P1')
    ax.plot([mmnt.dt for mmnt in mseq.vmmnt], [mmnt.p2 for mmnt in mseq.vmmnt], label='P2')
    ax.plot([mmnt.dt for mmnt in mseq.vmmnt], [mmnt.p3 for mmnt in mseq.vmmnt], label='P3')
    #ax.plot([mmnt.dt for mmnt in mseq.vmmnt], [mmnt.power for mmnt in mseq.vmmnt])
    #ax.legend(['P1', 'P2', 'P3', 'PWR'], loc='upper right')
    ax.legend(loc='upper right')
    ax.yaxis.set_minor_locator(mtckr.MultipleLocator(base=1))
    if mseq.td < timedelta(minutes=4): # Good for 1-4 min plots
      ax.xaxis.set_minor_locator(mdates.SecondLocator(bysecond=range(60)))
      ax.xaxis.set_major_locator(mdates.SecondLocator(bysecond=(0,15,30,45)))
      ax.xaxis.set_major_formatter(mdates.DateFormatter('%M:%S'))
    elif mseq.td < timedelta(minutes=16): # Good for the 10 min plots
      ax.xaxis.set_minor_locator(mdates.SecondLocator(bysecond=(10,20,30,40,50)))
      ax.xaxis.set_major_locator(mdates.MinuteLocator(byminute=range(60)))
      ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    elif mseq.td < timedelta(hours=16): # Good for the 10 h plots
      ax.xaxis.set_minor_locator(mdates.MinuteLocator(byminute=(10,20,30,40,50)))
      ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(24)))
      ax.xaxis.set_major_formatter(mdates.DateFormatter('%a %Hh'))
    elif mseq.td < timedelta(days=4): # Good for 1-4 day views
      ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=range(24)))
      ax.xaxis.set_major_locator(mdates.HourLocator(byhour=(0,6,12,18)))
      ax.xaxis.set_major_formatter(mdates.DateFormatter('%a %m-%d %H'))
    elif mseq.td < timedelta(days=16): # Good for 1-2 week views
      ax.xaxis.set_minor_locator(mdates.HourLocator(byhour=(6,12,18)))
      ax.xaxis.set_major_locator(mdates.WeekdayLocator(byweekday=range(7)))
      ax.xaxis.set_major_formatter(mdates.DateFormatter('%a %m-%d'))
    else: # Good for TD3600 month view (all data)
      ax.xaxis.set_minor_locator(mdates.WeekdayLocator(byweekday=range(7)))
      ax.xaxis.set_major_locator(mdates.WeekdayLocator(byweekday=(mdates.MO,mdates.WE,mdates.FR)))
      ax.xaxis.set_major_formatter(mdates.DateFormatter('%a %m-%d'))
    ax.xaxis.set_tick_params(which='major', length=6, width=4)
    ax.xaxis.set_tick_params(which='minor', length=4, width=2)
    ax.yaxis.set_tick_params(which='major', length=6, width=4)
    ax.yaxis.set_tick_params(which='minor', length=4, width=2)
    fig.autofmt_xdate()
    plt.savefig('out/%s.svg' % (mseq.pfdb.stem))
    fig.clear()
    plt.close(fig)

### HMM
elif action == 'hmm':
  ## Table printer function
  def print_table(tbl):
    for row in tbl: # {:>8}
      print(' '.join('{:>8}'.format(str(cell)) for cell in row))
  ## Implementation of the Viterbi algorithm:
  def viterbi(obs, states, start_p, trans_p, emit_p):
    V = [{}]
    for st in states:
      V[0][st] = {"prob": start_p[st] * emit_p[st][obs[0]], "prev": None}
    # Run Viterbi when t > 0
    for t in range(1, len(obs)):
      V.append({})
      for st in states:
        max_tr_prob = max(V[t-1][prev_st]["prob"]*trans_p[prev_st][st] for prev_st in states)
        for prev_st in states:
          if V[t-1][prev_st]["prob"] * trans_p[prev_st][st] == max_tr_prob:
            max_prob = max_tr_prob * emit_p[st][obs[t]]
            V[t][st] = {"prob": max_prob, "prev": prev_st}
            break
    # Print a table of steps from dictionary    
    tbl = [['steps'] + [str(i) for i in range(len(V))]]
    for state in V[0]:
      tbl.append([state] + ['{:.4f}'.format(v[state]["prob"]) for v in V])
    print('Computed table of steps:')
    print_table(tbl)
    opt = []
    # The highest probability
    max_prob = max(value["prob"] for value in V[-1].values())
    previous = None
    # Get most probable state and its backtrack
    for st, data in V[-1].items():
      if data["prob"] == max_prob:
        opt.append(st)
        previous = st
        break
    # Follow the backtrack till the first observation
    for t in range(len(V) - 2, -1, -1):
      opt.insert(0, V[t + 1][previous]["prev"])
      previous = V[t + 1][previous]["prev"]
    print('The most likely state sequence is:')
    print('  ' + ' -> '.join(opt))
    print('  Probability: {:.5f}'.format(max_prob))
  ## Data
  # Observations
  obs = ('jog', 'jog', 'clean', 'sleep', 'jog', 'sleep')
  # All possible (hidden) states
  states = ('sun', 'rain')
  # Initial start probabilities
  start_p = {'sun': 0.6, 'rain': 0.4}
  # Transition probabilities
  trans_p = {
    'sun' : {'sun': 0.7, 'rain': 0.3},
    'rain' : {'sun': 0.4, 'rain': 0.6}
  }
  # Emission probabilities
  emit_p = {
    'sun' : {'jog': 0.5, 'clean': 0.4, 'sleep': 0.1},
    'rain' : {'jog': 0.1, 'clean': 0.3, 'sleep': 0.6}
  }
  ## Execution
  # Run the algorithm
  viterbi(obs, states, start_p, trans_p, emit_p)

### HELP
else:
  print(__doc__)

sys.exit(0)
