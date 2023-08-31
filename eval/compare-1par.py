#!/usr/bin/env python3

import simplejson
import sys
import argparse

parser = argparse.ArgumentParser(
    description="Compare inferred with true parent phase")
parser.add_argument("-par", action="store", dest="par", type=str, default="dad")
parser.add_argument("-dir", action="store", dest="dir", type=str)

args = parser.parse_args()


def print_mismatch(parents, chr, marker_num, true, inf, codes):
  sys.stdout.write("Mismatch: {} chr{} marker {}: true {}|{}, inferred {}|{} ".
      format(parents, chr, marker_num, true[0], true[1], inf[0], inf[1]))
  #sys.stdout.write("          ")
  if codes != None:
    print(*codes, sep=', ')
  else:
    print("[None]")


if args.dir is None:
  args.dir = "2par-{}_miss".format(args.par)


orig = simplejson.loads(open("2par-orig/parent_phase.json", 'rb').read().decode('utf-8'))
inferred = simplejson.loads(open("{}/parent_phase.json".format(args.dir), 'rb').read().decode('utf-8'))

for parents, inf_data in inferred.items():
  if parents == "marker" or parents == "physpos" or parents == "chr" \
      or parents == "chrstr":
    continue

  par_idx = 0
  if args.par == "mom":
    par_idx = 1

  inf_phase = inf_data["parhaps"][par_idx]
  true_phase = orig[parents]["parhaps"][par_idx]

  # haplotype 0 in the inferred can correspond to haplotype 1 or 0
  # if this is 1, they're reversed, and to start, it's -1 for unknown
#  orientation = -1
#      if orientation < 0:
#        # No orientation set: try to assign
#        if inf[0] == inf[1] and inf[0] == '0':
#          # completely missing site: don't use to set orientation
#          continue
#
#        match = True
#        for hap in (0, 1):
#          if inf[hap] != true[0] and inf[hap] != '0':
#            match = False
#            break
#        if match:
#          orientation = 0
#          continue
#        # no match yet, try the other orientation:
#        for hap in (0, 1):      if orientation < 0:


  # TODO: comment
  # TODO: orientation is used inconsistently below: should be -1 for mismatch
  #       and None for uninformative.
  orientation = list()

  num_R_err = 0
  num_data_sites = [ 0, 0, 0 ]
  num_mismatches = 0

  # iterate over per-marker phases
  for idx, marker_data in enumerate(zip(inf_phase[0], inf_phase[1],
                                   inf_data["codes"], true_phase[0],
                                   true_phase[1], orig[parents]["codes"])):
    inf = marker_data[0:2]
    inf_codes = marker_data[2]
    true = marker_data[3:5]
    true_codes = marker_data[5]

    if inf_codes != None and inf_codes[0] == 'R':
      num_R_err += 1

    for codes in (true_codes, inf_codes):
      bad_site = False
      if codes != None and (codes[0] == '?' or codes[0] == 'E' \
          or codes[0] == 'R'):
        # bad phase in true or inferred data: skip
        bad_site = True
        orientation.append(None)
        break
    if bad_site:
      continue

    if inf[0] == '0' and inf[1] == '0':
      num_data_sites[0] += 1  # data for 0 alleles
    elif inf[0] == '0' or inf[1] == '0':
      num_data_sites[1] += 1  # data for 1 alleles
    else:
      num_data_sites[2] += 1  # data for both alleles

    if true[0] == true[1]:  # true is homozygous (or fully missing)
      if true[0] == '0':
        # true is missing: skip
        orientation.append(None)
        continue

      # true homozygous
      appended = False
      for hap in (0, 1):
        if inf[hap] != true[0] and inf[hap] != '0':
          chr = inferred["chr"][idx]
          print_mismatch(parents, chr, idx - inferred["chrstr"][chr], true, inf,
                         codes = inf_data["codes"][idx])
          num_mismatches += 1
          orientation.append(-1)
          appended = True
          break
      if not appended:
        orientation.append(None)
    elif true[0] == '0' or true[1] == '0': # true missing one allele
      # which true allele non-missing?
      non_miss_true = 0
      if true[0] == '0':
        non_miss_true = 1

      if inf[0] == inf[1]:
        # missing one allele and inferred is homozygous; will match in both
        # orientations or neither, so test:
        if inf[0] != true[non_miss_true] and inf[hap] != '0':
          chr = inferred["chr"][idx]
          print_mismatch(parents, chr, idx - inferred["chrstr"][chr], true, inf,
                         codes = inf_data["codes"][idx])
          num_mismatches += 1
          orientation.append(-1)
          appended = True
        else:
          orientation.append(None) # homozy inferred; half-true: uninformative
      elif inf[0] == '0' or true[1] == '0':
        # one missing site in both the true and inferred
        # which inferred allele non-missing?
        non_miss_inf = 0
        if inf[0] == '0':
          non_miss_inf = 1
        if inf[non_miss_inf] != true[non_miss_true]:
          # two non-missing alleles don't match, indicating
          # an orientation:
          # since the alleles mismatch, non_miss_inf ^ non_miss_true is the
          # opposite orientation, so:
          orientation.append(1 ^ non_miss_inf ^ non_miss_true)
          # The above can miss an inconsistent genotype, but not for biallelic
          # markers
        else:
          # the alleles match, but the site may be homozygous
          # so it's non-informative:
          orientation.append(None)
      else:
        # missing one allele and inferred is heterozygous
        # find orientation:
        if inf[non_miss_true] == true[non_miss_true]:
          orientation.append(0)
        elif inf[non_miss_true ^ 1] == true[non_miss_true]:
          orientation.append(1)
        else:
          # inconsistency:
          chr = inferred["chr"][idx]
          print_mismatch(parents, chr, idx - inferred["chrstr"][chr],
                         true, inf, codes = inf_data["codes"][idx])
          num_mismatches += 1
          orientation.append(-1)

    else:  # true is heterozygous
      if inf[0] == inf[1] and inf[0] == '0':
        # inferred is completely missing: no orientation
        orientation.append(None)
        continue

      match = True
      for hap in (0, 1):
        if inf[hap] != true[hap] and inf[hap] != '0':
          match = False
          break
      if match:
        orientation.append(0)
        continue
      # no match yet, try other orientation:
      match = True
      for hap in (0, 1):
        if inf[hap] != true[hap ^ 1] and inf[hap] != '0':
          match = False
          break
      if match:
        orientation.append(1)
        continue
      # inconsistency:
      chr = inferred["chr"][idx]
      print_mismatch(parents, chr, idx - inferred["chrstr"][chr],
                     true, inf, codes = inf_data["codes"][idx])
      num_mismatches += 1
      orientation.append(-1)

  if len(orientation) != len(inferred["chr"]):
    print("ERROR: length of orientation list different from number of 'chr' entries? {} and {}".format(len(orientation), len(inferred["chr"])))
    sys.exit(1)

  last_chr = inferred["chr"][0]
  last_orient = None
  last_orient_idx = 0
  num_switches = 0
  the_switches = list()
  for idx, (orient, chr) in enumerate(zip(orientation, inferred["chr"])):
    if last_chr != chr:
      # end of chromosome, print number of switches and reset
      print("{}, chr {}: {} {}".format(parents, last_chr, num_switches, ",".join(str(x) for x in the_switches)))
      last_chr = chr
      last_orient = None
      last_orient_idx = idx
      num_switches = 0
      the_switches = list()
        
    if orient != None and orient >= 0:
      if orient != last_orient:
        if last_orient != None:
          num_switches += 1
          the_switches.append(idx - last_orient_idx)
          idx = last_orient_idx
        last_orient = orient

  print("{}, chr {}: {} {}".format(parents, last_chr, num_switches, ",".join(str(x) for x in the_switches)))
  print("STAT R sites    {}: {}".format(parents, num_R_err))
  total_sites = sum(num_data_sites)
  print("STAT Mismatches {}: {:.6f}".format(parents, num_mismatches / total_sites))
  print("STAT Full miss  {}: {:.6f}".format(parents, num_data_sites[0] / total_sites))
  print("STAT Half data  {}: {:.6f}".format(parents, num_data_sites[1] / total_sites))
  print("STAT Full data  {}: {:.6f}".format(parents, num_data_sites[2] / total_sites))


