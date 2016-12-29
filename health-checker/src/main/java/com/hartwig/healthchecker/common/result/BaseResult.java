package com.hartwig.healthchecker.common.result;

import java.io.Serializable;

import com.hartwig.healthchecker.common.checks.CheckType;

import org.jetbrains.annotations.NotNull;

public interface BaseResult extends Serializable {

    long serialVersionUID = -4752339157661751000L;

    @NotNull
    CheckType getCheckType();
}
