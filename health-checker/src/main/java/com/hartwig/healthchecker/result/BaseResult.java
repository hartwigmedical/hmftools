package com.hartwig.healthchecker.result;

import java.io.Serializable;

import com.hartwig.healthchecker.runners.CheckType;

import org.jetbrains.annotations.NotNull;

public interface BaseResult extends Serializable {

    long serialVersionUID = -4752339157661751000L;

    @NotNull
    CheckType getCheckType();
}
