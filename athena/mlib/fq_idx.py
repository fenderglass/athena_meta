import os
import sys
import cPickle as pickle
from itertools import groupby

from Bio import bgzf

import util

class FastqIndex(object):

  file_suffix = '.fqidx.p'

  @staticmethod
  def get_index_path(fq_path):
    return fq_path + FastqIndex.file_suffix

  @property
  def bcode_set(self):
    if self._bcode_set == None:
      self._bcode_set = set(self._bcode_off_map.keys())
    return self._bcode_set

  def __init__(
    self,
    fq_path,
    logger=None,
  ):
    
    self.logger = logger
    self.fq_path = fq_path
    self.index_path = self.get_index_path(fq_path)
    self._bcode_set = None
    self._bcode_off_map = None

    if not os.path.isfile(self.index_path):
      self.__build_index__()
    else:
      self.__load_index__()

    self.f_map = None
    self.open()

  def open(self):
    assert self.f_map == None, "fp map already populated"
    self.f_map = {}

    if self.fq_path.endswith('.gz'):
      index_name = self.fq_path + "i"
      if not os.path.exists(index_name):
        raise Exception("Only BGZF compression is supported")

      handle = bgzf.BgzfReader(self.fq_path)
      self.gzipped = True
    else:
      handle = open(self.fq_path)
      self.gzipped = False

    self.f_map[self.fq_path] = handle
    return self

  def close(self):
    for f in self.f_map.values():
      f.close()
    return

  def __enter__(self):
    return self

  def __exit__(self, exc_type, exc_value, traceback):
    self.close()

  def __build_index__(self):  
    numbytes = 0
    self._bcode_off_map = {}
    num_pe = 0

    if self.fq_path.endswith('.gz'):
      index_name = self.fq_path + "i"
      if not os.path.exists(index_name):
        raise Exception("Only BGZF compression is supported")

      handle = bgzf.BgzfReader(self.fq_path)
    else:
      handle = open(self.fq_path)

    seen_set = set()
    for bcode, reads_iter in groupby(
      util.fastq_iter_pos(handle),
      lambda(x): x[0],
    ):
      assert bcode == None or bcode not in seen_set, \
"fastq {} NOT in barcode sorted order. Ensure reads that share barcodes \
are in a block together".format(self.fq_path)
      seen_set.add(bcode)
      for _, qname, file_pos, lines in reads_iter:
        if bcode != None and bcode not in self._bcode_off_map:
          self._bcode_off_map[bcode] = file_pos
        num_pe += 1
    handle.close()

    num_bcodes = len(filter(
      lambda(b): b.endswith('-1'),
      self._bcode_off_map.keys(),
    ))
    assert num_bcodes > 0, \
      "no barcodes specified in fastq {}".format(self.fq_path)
    self.logger.log('fqinfo${},{},{}'.format(
      num_pe, len(self._bcode_off_map), num_bcodes,
    ))
    print 'writing index for fqs'
    for fq_path in [self.fq_path]:
      print '  -', fq_path
    util.write_pickle(self.index_path, self._bcode_off_map)

  def __load_index__(self):  
    self._bcode_off_map = util.load_pickle(self.index_path)

  def get_reads(self, bcode):
    assert self._bcode_off_map != None, 'index {} not loaded'.format(self.index_path)

    if bcode not in self.bcode_set:
      raise StopIteration

    offset = self._bcode_off_map[bcode]
    f = self.f_map[self.fq_path]
    if not self.gzipped:
      f.seek(offset, 0)
    else:
      f.seek(offset)

    for e in util.fastq_iter(f):
      _bcode = e[0]
      if _bcode != bcode:
        break
      yield e

    raise StopIteration

