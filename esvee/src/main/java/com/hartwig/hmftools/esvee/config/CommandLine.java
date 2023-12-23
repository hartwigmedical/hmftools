package com.hartwig.hmftools.esvee.config;

import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;

@Retention(RetentionPolicy.RUNTIME)
public @interface CommandLine {
    /** If not supplied, defaults to the field name converted to snake_case. */
    String name() default "";
    String description();
}
