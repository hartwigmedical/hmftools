package com.hartwig.hmftools.esvee.config;

import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;

@Retention(RetentionPolicy.RUNTIME)
public @interface Optional {
    String defaultValue();
}
