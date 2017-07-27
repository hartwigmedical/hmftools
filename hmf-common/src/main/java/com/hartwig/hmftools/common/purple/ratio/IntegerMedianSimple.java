package com.hartwig.hmftools.common.purple.ratio;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;

class IntegerMedianSimple {

    private int count;
    private List<Integer> list = Lists.newArrayList();

    public void addRead(int read) {
        count++;
        list.add(read);
    }

    public int median() {
        if (count > 0) {
            Collections.sort(list);
            return list.size() % 2 == 0 ? (list.get(count / 2) + list.get(count / 2 - 1)) / 2 : list.get(count / 2);
        }
        return 0;
    }

}
