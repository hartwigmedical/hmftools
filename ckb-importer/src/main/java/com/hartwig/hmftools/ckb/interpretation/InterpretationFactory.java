package com.hartwig.hmftools.ckb.interpretation;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.datamodel.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.datamodelinterpretation.CkbEntryInterpretation;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class InterpretationFactory {

    private InterpretationFactory() {

    }

    private static final Logger LOGGER = LogManager.getLogger(InterpretationFactory.class);

    public static List<CkbEntryInterpretation> interpretationCkbDataModel(@NotNull CkbJsonDatabase ckbEntry) {
        List<CkbEntryInterpretation> CkbEntryInterpretation = Lists.newArrayList();
        return CkbEntryInterpretation;
    }

}
