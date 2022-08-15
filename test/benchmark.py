# -*- coding: utf-8 -*-

import timeit
import numpy as np
import pandas as pd
from pymatch2.MatcherPlus import MatcherPlus

class MatcherCmp:
    def __init__(self, col_identity = "identity", col_treatment = "treatment", col_matchid = "match_id",
            col_scores = "scores"):
        self.col_identity_ = col_identity
        self.col_treatment_ = col_treatment
        self.col_matchid_ = col_matchid
        self.col_scores_ = col_scores

    def _merged(self, matched):
        reduced = matched[[self.col_matchid_, self.col_identity_, self.col_treatment_, self.col_scores_]]
        mask = (reduced[self.col_treatment_] == True)
        treated = reduced[mask]
        control = reduced[~mask]
        merged = pd.merge(treated, control, how = "inner", on = self.col_matchid_, suffixes = ("_x", "_y"))
        col_diff = "diff"
        merged[col_diff] = abs(merged[self.col_scores_ + "_x"] - merged[self.col_scores_ + "_y"])
        merged[self.col_identity_] = merged[self.col_identity_ + "_x"]
        return merged[[self.col_identity_, col_diff]], col_diff

    def _compare(self, lhs, rhs):
        assert lhs.shape == rhs.shape
        assert lhs.shape[0] > 0
        lmerged, col_diff = self._merged(lhs)
        rmerged, _ = self._merged(rhs)
        merged = pd.merge(lmerged, rmerged, how = "inner", on = self.col_identity_, suffixes = ("_x", "_y"))
        col_did = "did"
        merged[col_did] = merged[col_diff + "_x"] - merged[col_diff + "_y"]
        return merged, col_did

    def less_than(self, lhs, rhs):
        merged, col_did = self._compare(lhs, rhs)
        return np.all(merged[col_did] < 0)

    def lt_pct(self, lhs, rhs):
        merged, col_did = self._compare(lhs, rhs)
        return np.mean(merged[col_did] < 0), merged.shape

    def less_eql(self, lhs, rhs):
        merged, col_did = self._compare(lhs, rhs)
        return np.all(merged[col_did] <= 0)

    def equal(self, lhs, rhs):
        merged, col_did = self._compare(lhs, rhs)
        return np.all(merged[col_did] == 0)

class MatcherTest:
    def __init__(self, ntest = 10000, nctrl = 100000):
        self.ntest_ = ntest
        self.nctrl_ = nctrl
        self.npopu_ = ntest + nctrl

        self.col_treatment_ = "treatment"
        col_scores = "scores"
        self.df_ = pd.DataFrame({
            self.col_treatment_: [0] * self.npopu_,
            col_scores: np.random.random(self.npopu_)
        })

        itest = np.random.choice(self.npopu_, ntest, replace = False)
        self.df_.loc[itest, self.col_treatment_] = 1

        self.col_identity_ = "identity"
        self.df_[self.col_identity_] = self.df_.index.values

        test_mask = (self.df_[self.col_treatment_] == 1)
        self.matcher_ = MatcherPlus(self.df_[test_mask], self.df_[~test_mask], self.col_treatment_)

        col_matchid = "match_id"
        self.compare_ = MatcherCmp(self.col_identity_, self.col_treatment_, col_matchid, col_scores)

    def timed_match(self, method = "nnm", nmatches = 1, threshold = 0.001):
        self.matcher_.data = self.df_.copy()
        start = timeit.default_timer()
        if method == "nnm":
            self.matcher_.match_nnm(nmatches = nmatches)
        elif method == "bs":
            self.matcher_.match_bs(nmatches = nmatches)
        elif method == "merge":
            self.matcher_.match_merge(nmatches = nmatches)
        else:
            self.matcher_.match(method = method, nmatches = nmatches, threshold = threshold)
        end = timeit.default_timer()
        return end - start, self.matcher_.matched_data

    def assert_lt(self, nmatches = 1, threshold = 0.001):
        """ nnm match less than random match """
        rnd_time, rnd_matched = self.timed_match("random", nmatches = nmatches, threshold = threshold)
        nnm_time, nnm_matched = self.timed_match("nnm", nmatches = nmatches, threshold = threshold)
        assert self.compare_.less_than(nnm_matched, rnd_matched)
        ## assert nnm_time <= rnd_time

        print(self.ntest_, self.nctrl_, rnd_time, "-", nnm_time, "-")

    def assert_le(self, nmatches = 1, threshold = 0.001):
        """ nnm match less than or equal to random match """
        rnd_time, rnd_matched = self.timed_match("random", nmatches = nmatches, threshold = threshold)
        nnm_time, nnm_matched = self.timed_match("nnm", nmatches = nmatches, threshold = threshold)
        assert self.compare_.less_eql(nnm_matched, rnd_matched)
        ## assert nnm_time <= rnd_time

        bs_time, bs_matched = self.timed_match("bs", nmatches = nmatches, threshold = threshold)
        assert self.compare_.less_eql(bs_matched, rnd_matched)

        mt_time, mt_matched = self.timed_match("merge", nmatches = nmatches, threshold = threshold)
        assert self.compare_.less_eql(mt_matched, rnd_matched)

        print(self.ntest_, self.nctrl_, rnd_time, "-", nnm_time, bs_time, mt_time)

    def assert_eq(self, nmatches = 1, threshold = 0.001):
        """ nnm match equal to min match """
        min_time, min_matched = self.timed_match("min", nmatches = nmatches, threshold = threshold)
        nnm_time, nnm_matched = self.timed_match("nnm", nmatches = nmatches, threshold = threshold)
        assert self.compare_.equal(nnm_matched, min_matched)
        # assert nnm_time <= min_time

        print(self.ntest_, self.nctrl_, "-", min_time, nnm_time, "-")

    def assert_match(self, nmatches = 1, threshold = 0.001):
        rnd_time, rnd_matched = self.timed_match("random", nmatches = nmatches, threshold = threshold)
        nnm_time, nnm_matched = self.timed_match("nnm", nmatches = nmatches, threshold = threshold)
        assert self.compare_.less_eql(nnm_matched, rnd_matched)

        min_time, min_matched = self.timed_match("min", nmatches = nmatches, threshold = threshold)
        assert self.compare_.equal(nnm_matched, min_matched)

        bs_time, bs_matched = self.timed_match("bs", nmatches = nmatches, threshold = threshold)
        assert self.compare_.equal(nnm_matched, bs_matched)

        mt_time, mt_matched = self.timed_match("merge", nmatches = nmatches, threshold = threshold)
        assert self.compare_.equal(nnm_matched, mt_matched)

        print(self.ntest_, self.nctrl_, rnd_time, min_time, nnm_time, bs_time, mt_time)

if __name__ == "__main__":
    np.random.seed(1024)
    print("# test", "# control", "random", "min", "nnm", "bs", "merge", flush = True)
    for ntest, nctrl in [(10000, 200000), (15000, 200000), (20000, 200000)]:
        test = MatcherTest(ntest, nctrl)
        test.assert_match()

    for ntest, nctrl in [(50000, 900000), (100000, 1000000), (500000, 5000000)]:
        test = MatcherTest(ntest, nctrl)
        test.assert_le()
