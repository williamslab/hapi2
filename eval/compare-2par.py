#!/usr/bin/env python3

import simplejson
import sys
import argparse

parser = argparse.ArgumentParser(
    description="Compare inferred with true parent genotypes and phase")
parser.add_argument("-true_dir", action="store", dest="true_dir", default='2par-orig', type=str)
parser.add_argument("-inf_dir", action="store", dest="inf_dir", default='2par-both_miss', type=str)
parser.add_argument("-parent", action="store", dest="parent", default=2, type=int)

args = parser.parse_args()

# how many parent orientations? 2 if we're analyzing both parents, 1 otherwise
# this is for P regions -- which parent is the dad and which is the mom can
# flip, so we examine both possibilities
numParentOrient = 2
if args.parent != 2:
  numParentOrient = 1

# which parents are we comparing?
parentsToExamine = [0, 1]
if args.parent != 2:
  parentsToExamine = [args.parent]


def finishPregionOrChrom(parents, chr, mismatchRecords, the_orientation,
                         phase_orientation):
  betterParOrient = 0
  if numParentOrient == 2 and len(mismatchRecords[1]) < len(mismatchRecords[0]):
    betterParOrient = 1
  flushMismatches(parents, chr, mismatchRecords[betterParOrient])
  for parent in parentsToExamine:
    the_orientation[parent].extend(
              phase_orientation[betterParOrient][parent] )
  for po in range(numParentOrient):
    for parent in parentsToExamine:
      phase_orientation[po][parent].clear()

  lowerTally = len(mismatchRecords[betterParOrient])
  mismatchRecords[0].clear()
  mismatchRecords[1].clear()

  return lowerTally


def flushMismatches(parents, chr, bestMismatchRecords):
  for r in bestMismatchRecords:
    print_mismatch(parents, chr, r[0], r[1], r[2], r[3], r[4])

def print_mismatch(parents, chr, trueParent, marker_num, true, inf, codes):
  sys.stdout.write(
      "Mismatch: {}-chr{} parent {} marker {}: true {}|{}, inf {}|{} ".
      format(parents, chr, trueParent, marker_num, true[0], true[1], inf[0],
             inf[1]))
  #sys.stdout.write("          ")
  if codes != None:
    print(*codes, sep=', ')
  else:
    print("[None]")

if args.parent < 0 or args.parent > 2:
  print("Error: -parent can be 0, 1, or 2 for both ({} given)".format(args.parent))
  sys.exit(2)

orig = simplejson.loads(open("{}/parent_phase.json".format(args.true_dir), 'rb').read().decode('utf-8'))
inferred = simplejson.loads(open("{}/parent_phase.json".format(args.inf_dir), 'rb').read().decode('utf-8'))

for parentsStr, inf_data in inferred.items():
  if parentsStr == "marker" or parentsStr == "physpos" or parentsStr == "chr"\
      or parentsStr == "chrstr":
    continue

  # dad and mom's phase
  true_phase = orig[parentsStr]["parhaps"]
  inf_phase =  inf_data["parhaps"]

  # TODO: orientation is used inconsistently below: should be -1 for mismatch
  #       and None for uninformative.
  # we store the phase orientation relative to the truth (see below)
  # when both parents are missing, we don't know which parent is
  # which, so <phase_orientation> stores two entries corresponding
  # to the different orientations; we append the preferred orientation
  # to <the_orientation> once we know which one is
  phase_orientation = [ [ list(), list() ], [ list(), list() ] ]
  the_orientation = [ list(), list() ]

  num_R_err = 0
  # TODO: actually, there are "parent" switch errors here; ideally would do
  # something like that done with orientation
  num_data_sites = [ [ 0, 0, 0 ], [ 0, 0, 0 ] ]

  mismatchRecords = [ list(), list() ]
  bestMismatchTally = 0

  prevChr = -1
  mostRecentP = -1

  marker_idx = -1
  while marker_idx < len(inf_data["codes"]) - 1:
    marker_idx += 1

    curChr = inferred["chr"][marker_idx]
    inf_codes = inf_data["codes"][marker_idx]
    true_codes = orig[parentsStr]["codes"][marker_idx]

    if curChr != prevChr:
      if prevChr != -1:
        ret = finishPregionOrChrom(parentsStr, prevChr, mismatchRecords,
                                   the_orientation, phase_orientation)
        bestMismatchTally += ret
      mostRecentP = -1
      prevChr = curChr
      # TODO: maybe rm after more testing
      assert(phase_orientation == [ [ list(), list() ], [ list(), list() ] ])
 
    if inf_codes != None and (inf_codes[0] == 'P' or inf_codes[0] == 'PC' or \
        inf_codes[0] == 'PA'):
      if mostRecentP >= 0:
        ret = finishPregionOrChrom(parentsStr, prevChr, mismatchRecords,
                                   the_orientation, phase_orientation)
        bestMismatchTally += ret
      mostRecentP = marker_idx

    if inf_codes != None and inf_codes[0] == 'R':
      num_R_err += 1
      if true_codes == None or true_codes[0] != 'R':
        print("Unnecessary R: {}-chr{} marker {}".
                         format(parentsStr, curChr,
                                marker_idx - inferred["chrstr"][curChr]))

    bad_site = False
    for codes in (true_codes, inf_codes):
      if codes != None and (codes[0] == '?' or codes[0] == 'E' \
          or codes[0] == 'R'):
        # bad phase in true or inferred data: skip
        bad_site = True
        for po in range(numParentOrient):
          phase_orientation[po][0].append(None)
          phase_orientation[po][1].append(None)
        break
    if bad_site:
      continue

    true = [ [ true_phase[0][0][marker_idx], true_phase[0][1][marker_idx] ],
             [ true_phase[1][0][marker_idx], true_phase[1][1][marker_idx] ] ]

    # for numParentOrient == 2, we don't know which parent is which, so we'll
    # try both orientations
    for par_orient in range(numParentOrient):
      inf = [ [ inf_phase[0 ^ par_orient][0][marker_idx],
                inf_phase[0 ^ par_orient][1][marker_idx] ],
              [ inf_phase[1 ^ par_orient][0][marker_idx],
                inf_phase[1 ^ par_orient][1][marker_idx] ] ]

      for parent in parentsToExamine: #compare inferred and true for each parent

        # count number of sites of various kinds
        # keep only for the "default" orientation (can swap at R markers if
        # needed)
        if par_orient == 0:
          if inf[parent][0] == '0' and inf[parent][1] == '0':
            num_data_sites[parent][0] += 1  # data for 0 alleles
          elif inf[parent][0] == '0' or inf[parent][1] == '0':
            num_data_sites[parent][1] += 1  # data for 1 alleles
          else:
            num_data_sites[parent][2] += 1  # data for both alleles

        if true[parent][0] == true[parent][1]:
          # true is homozygous (or fully missing)
          if true[parent][0] == '0':
            # true is missing: skip
            phase_orientation[par_orient][parent].append(None)
            continue

          # true homozygous
          appended = False
          for hap in (0, 1):
            if inf[parent][hap] != true[parent][0] and inf[parent][hap] !='0':
              mismatchRecords[par_orient].append(
                  (parent, marker_idx - inferred["chrstr"][curChr],
                   true[parent], inf[parent], inf_codes) )
              phase_orientation[par_orient][parent].append(-1)
              appended = True
              break
          if not appended:
            phase_orientation[par_orient][parent].append(None)
        elif true[parent][0] == '0' or true[parent][1] == '0':
          # true missing one allele
          # which true allele non-missing?
          non_miss_true = 0
          if true[parent][0] == '0':
            non_miss_true = 1

          if inf[parent][0] == inf[parent][1]:
            # missing one allele and inferred is homozygous; will match in
            # both orientations or neither, so test:
            if inf[parent][0] != true[parent][non_miss_true] and inf[parent][hap] != '0':
              mismatchRecords[par_orient].append(
                  (parent, marker_idx - inferred["chrstr"][curChr],
                   true[parent], inf[parent], inf_codes) )
              phase_orientation[par_orient][parent].append(-1)
              appended = True
            else:
              # homozy inferred; half-true: uninformative
              phase_orientation[par_orient][parent].append(None)
          elif inf[parent][0] == '0' or true[parent][1] == '0':
            # one missing site in both the true and inferred
            # which inferred allele non-missing?
            non_miss_inf = 0
            if inf[parent][0] == '0':
              non_miss_inf = 1
            if inf[parent][non_miss_inf] != true[parent][non_miss_true]:
              # two non-missing alleles don't match, indicating
              # an orientation:
              # since the alleles mismatch, non_miss_inf ^ non_miss_true is
              # the opposite orientation, so:
              phase_orientation[par_orient][parent].append(
                  1 ^ non_miss_inf ^ non_miss_true)
              # The above can miss an inconsistent genotype, but not for
              # biallelic markers
            else:
              # the alleles match, but the site may be homozygous
              # so it's non-informative:
              phase_orientation[par_orient][parent].append(None)
          else:
            # missing one allele and inferred is heterozygous
            # find orientation:
            if inf[parent][non_miss_true] == true[parent][non_miss_true]:
              phase_orientation[par_orient][parent].append(0)
            elif inf[parent][non_miss_true ^ 1] == true[parent][non_miss_true]:
              phase_orientation[par_orient][parent].append(1)
            else:
              # inconsistency:
              mismatchRecords[par_orient].append(
                  (parent, marker_idx - inferred["chrstr"][curChr],
                   true[parent], inf[parent], inf_codes) )
              phase_orientation[par_orient][parent].append(-1)

        else:  # true is heterozygous
          if inf[parent][0] == inf[parent][1] and inf[parent][0] == '0':
            # inferred is completely missing: no orientation
            phase_orientation[par_orient][parent].append(None)
            continue

          match = True
          for hap in (0, 1):
            if inf[parent][hap] != true[parent][hap] and \
                inf[parent][hap] != '0':
              match = False
              break
          if match:
            phase_orientation[par_orient][parent].append(0)
            continue
          # no match yet, try other orientation:
          match = True
          for hap in (0, 1):
            if inf[parent][hap] != true[parent][hap ^ 1] and \
                inf[parent][hap] != '0':
              match = False
              break
          if match:
            phase_orientation[par_orient][parent].append(1)
            continue
          # inconsistency:
          mismatchRecords[par_orient].append(
              (parent, marker_idx - inferred["chrstr"][curChr],
               true[parent], inf[parent], inf_codes) )
          phase_orientation[par_orient][parent].append(-1)

  # finished analyzing all chromosomes: put last set of mismatches and
  # orientations into their permanent containers, considering all parent
  # orientations
  ret = finishPregionOrChrom(parentsStr, prevChr, mismatchRecords,
                             the_orientation, phase_orientation)
  bestMismatchTally += ret

  for parent in parentsToExamine:
    if len(the_orientation[parent]) != len(inferred["chr"]):
      print("ERROR: length of orientation list different from number of 'chr' entries? {} and {}".format(len(the_orientation[parent]), len(inferred["chr"])))
      sys.exit(1)

  last_chr = inferred["chr"][0]
  last_orient = [None, None]
  last_orient_idx = 0
  num_switches = 0
  the_switches = list()
  for idx, (p0orient, p1orient, chr) in enumerate(zip(the_orientation[0],
                                                      the_orientation[1],
                                                      inferred["chr"])):
    orient = (p0orient, p1orient)

    if last_chr != chr:
      # end of chromosome, print number of switches and reset
      print("{}, chr {}: {} {}".format(parentsStr, last_chr, num_switches, ",".join(str(x) for x in the_switches)))
      last_chr = chr
      last_orient = [None, None]
      last_orient_idx = idx
      num_switches = 0
      the_switches = list()
        
    for parent in parentsToExamine:
      if orient[parent] != None and orient[parent] >= 0:
        if orient[parent] != last_orient[parent]:
          if last_orient[parent] != None:
            num_switches += 1
            the_switches.append(idx - last_orient_idx)
            idx = last_orient_idx
          last_orient[parent] = orient[parent]

  print("{}, chr {}: {} {}".format(parentsStr, last_chr, num_switches, ",".join(str(x) for x in the_switches)))
  print("STAT R sites    {}: {}".format(parentsStr, num_R_err))
  total_sites = []
  mismatchDenom = 0
  for parent in parentsToExamine:
    parTot = sum(num_data_sites[parent])
    total_sites.append( parTot )
    mismatchDenom += parTot
  if len(parentsToExamine) == 1:
    num_data_sites[0] = num_data_sites[ parentsToExamine[0] ]
  print("STAT Mismatches {}: {:.6f}".format(parentsStr,
    bestMismatchTally / mismatchDenom))
  print("STAT Missing    {}: {:.6f}".format(parentsStr,
        num_data_sites[0][0] / total_sites[0]),
        end='')
  if len(parentsToExamine) > 1:
    print(", {:.6f}".format(num_data_sites[1][0] / total_sites[1]))
  else:
    print()

  print("STAT Half data  {}: {:.6f}".format(parentsStr,
        num_data_sites[0][1] / total_sites[0]),
        end='')
  if len(parentsToExamine) > 1:
    print(", {:.6f}".format(num_data_sites[1][1] / total_sites[1]))
  else:
    print()

  print("STAT Full data  {}: {:.6f}".format(parentsStr,
        num_data_sites[0][2] / total_sites[0]),
        end='')
  if len(parentsToExamine) > 1:
    print(", {:.6f}".format(num_data_sites[1][2] / total_sites[1]))
  else:
    print()


