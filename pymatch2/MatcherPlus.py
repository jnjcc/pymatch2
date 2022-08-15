# -*- coding: utf-8 -*-

import pandas as pd
from pymatch.Matcher import Matcher

class MatcherPlus(Matcher):
    """ Matcher Class with faster Nearest Neighbor Matching algorithms """

    def __init__(self, test, control, yvar, formula = None, exclude = []):
        super().__init__(test, control, yvar, formula, exclude)
        self.nmatches_ = 1
        self.indices_ = []
        self.match_ids_ = []
        self.cur_match_id_ = 0

    def match_nnm(self, nmatches = 1):
        """ NOTICE: nmatches = min{`nmatches', (l2 - l1) + (r2 - r1)} """
        if "scores" not in self.data.columns:
            print("Propensity Scores have not been calculated. Using defaults...")
            self.fit_scores()
            self.predict_scores()

        df_data = self.data.sort_values("scores").reset_index(drop = True)
        self.indices_ = [] ## indices of `df_data' after matching
        self.match_ids_ = []
        self.cur_match_id_ = 0
        self.nmatches_ = nmatches

        EOS = df_data.shape[0]
        state = 0
        ## (l1, l2], (l2, r1], (r1, r2]
        l1, l2, r1, r2 = [-1] * 4
        for i in range(EOS + 1):
            ## 1 if in treated group, 0 if in control group, -1 if EOS
            isymbol = int(df_data.iloc[i][self.yvar]) if i < EOS else -1
            if state == 0:
                if isymbol == 0:
                    l2 = i
                    state = 1
                elif isymbol == 1:
                    r1 = i
                    state = 2
            elif state == 1:
                if isymbol == 0:
                    l2 = i
                    state = 1
                elif isymbol == 1:
                    r1 = i
                    state = 2
            elif state == 2:
                if isymbol == 1:
                    r1 = i
                    state = 2
                elif isymbol == 0:
                    r2 = i
                    state = 3
                elif isymbol == -1:
                    self._match_interval(df_data, l1, l2, r1, r2)
            elif state == 3:
                if isymbol == 0:
                    r2 = i
                    state = 3
                elif isymbol == 1:
                    self._match_interval(df_data, l1, l2, r1, r2)
                    l1 = r1
                    l2 = r2
                    r1 = i
                    state = 2
                elif isymbol == -1:
                    self._match_interval(df_data, l1, l2, r1, r2)
        self.matched_data = df_data.loc[self.indices_]
        self.matched_data["match_id"] = self.match_ids_
        self.matched_data["record_id"] = self.matched_data.index

    def _match_interval(self, df_data, l1, l2, r1, r2):
        itreat = l2 + 1
        while itreat <= r1:
            iscore = df_data.iloc[itreat]["scores"]
            indices = self._merge(iscore, df_data, l1, l2, r1, r2)
            indices.append(df_data.index[itreat])
            self.indices_.extend(indices)
            self.match_ids_.extend([self.cur_match_id_] * len(indices))
            self.cur_match_id_ += 1
            itreat += 1

    def _merge(self, iscore, dframe, l1, l2, r1, r2):
        """ (l1, l2], (r1, r2] """
        indices = []
        lidx = l2
        ridx = r1 + 1
        while (len(indices) < self.nmatches_) and ((lidx > l1) and (ridx <= r2)):
            ldiff = iscore - dframe.iloc[lidx]["scores"]
            rdiff = dframe.iloc[ridx]["scores"] - iscore
            if ldiff < rdiff:
                indices.append(dframe.index[lidx])
                lidx -= 1
            elif ldiff >= rdiff:
                indices.append(dframe.index[ridx])
                ridx += 1
        while (len(indices) < self.nmatches_) and (lidx > l1):
            indices.append(dframe.index[lidx])
            lidx -= 1
        while (len(indices) < self.nmatches_) and (ridx <= r2):
            indices.append(dframe.index[ridx])
            ridx += 1
        return indices

    def match_bs(self, nmatches = 1):
        """ Matching using binary search algorithm """
        if "scores" not in self.data.columns:
            print("Propensity Scores have not been calculated. Using defaults...")
            self.fit_scores()
            self.predict_scores()

        df_test = self.data[self.data[self.yvar] == True]
        df_ctrl = self.data[self.data[self.yvar] == False]
        df_ctrl = df_ctrl.sort_values("scores").reset_index(drop = True)
        self.indices_ = [] ## indices of `df_ctrl' after matching
        self.match_ids_ = list(range(df_test.shape[0])) ## match_ids of `matched_data'
        self.nmatches_ = nmatches

        ntest = df_test.shape[0]
        nctrl = df_ctrl.shape[0]
        for i in range(ntest):
            iscore = df_test.iloc[i]["scores"]
            match = self._binary_search(iscore, df_ctrl, 0, nctrl)
            eidx = (df_ctrl.shape[0] - 1)
            indices = [df_ctrl.index[match]]
            if self.nmatches_ > 1:
                ## it must be true that (scores[match-1] <= iscore <= scores[match+1])
                ##   thus (-1, match-1], (match, eidx]
                merged = self._merge(iscore, df_ctrl, -1, match - 1, match, eidx)
                indices.extend(merged)
            self.indices_.extend(indices)
            self.match_ids_.extend([i] * len(indices))
        self.matched_data = pd.concat([df_test, df_ctrl.loc[self.indices_]], axis = 0)
        self.matched_data["match_id"] = self.match_ids_
        self.matched_data["record_id"] = self.matched_data.index

    def _binary_search(self, iscore, df_ctrl, left, right):
        """ Binary search in `df_ctrl' for Nearest Neighbor of `iscore' """
        l = left
        r = right - 1
        match = left
        mdiff = float("inf")
        while l <= r:
            mid = l + (r - l) // 2
            diff = df_ctrl.iloc[mid]["scores"] - iscore
            adiff = abs(diff)
            if adiff < mdiff:
                mdiff = adiff
                match = mid
            if diff > 0:
                r = mid - 1
            elif diff < 0:
                l = mid + 1
            else:
                match = mid
                l = mid
                r = mid
                break
        return match

    def match_merge(self, nmatches = 1):
        """ Matching using merge-sort alike algorithm """
        if "scores" not in self.data.columns:
            print("Propensity Scores have not been calculated. Using defaults...")
            self.fit_scores()
            self.predict_scores()

        df_test = self.data[self.data[self.yvar] == True]
        df_ctrl = self.data[self.data[self.yvar] == False]
        df_test = df_test.sort_values("scores").reset_index(drop = True)
        df_ctrl = df_ctrl.sort_values("scores").reset_index(drop = True)
        self.indices_ = [] ## indices of `df_ctrl' after matching
        self.match_ids_ = list(range(df_test.shape[0])) ## match_ids of `matched_data'
        self.nmatches_ = nmatches

        ntest = df_test.shape[0]
        nctrl = df_ctrl.shape[0]
        jctrl = 0
        for i in range(ntest):
            iscore = df_test.iloc[i]["scores"]
            while (jctrl < nctrl) and (df_ctrl.iloc[jctrl]["scores"] < iscore):
                jctrl += 1
            eidx = (df_ctrl.shape[0] - 1)
            ## it must be true that (scores[jctrl-1] <= iscore <= scores[jctrl])
            ##   thus (-1, jctrl-1], (jctrl-1, eidx]
            indices = self._merge(iscore, df_ctrl, -1, jctrl - 1, jctrl - 1, eidx)
            self.indices_.extend(indices)
            self.match_ids_.extend([i] * len(indices))

        self.matched_data = pd.concat([df_test, df_ctrl.loc[self.indices_]], axis = 0)
        self.matched_data["match_id"] = self.match_ids_
        self.matched_data["record_id"] = self.matched_data.index
