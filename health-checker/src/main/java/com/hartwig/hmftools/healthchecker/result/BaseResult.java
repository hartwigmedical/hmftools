package com.hartwig.hmftools.healthchecker.result;

import java.io.Serializable;

import com.hartwig.hmftools.healthchecker.runners.CheckType;

import org.jetbrains.annotations.NotNull;

public interface BaseResult extends Serializable {

    long serialVersionUID = -4752339157661751000L;

    @NotNull
    CheckType getCheckType();
}
