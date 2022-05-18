package com.hartwig.hmftools.orange.report.util;

import java.net.MalformedURLException;

import com.itextpdf.io.IOException;
import com.itextpdf.io.image.ImageDataFactory;
import com.itextpdf.layout.element.Image;

import org.jetbrains.annotations.NotNull;

public final class Images {

    private Images() {
    }

    @NotNull
    public static Image build(@NotNull String path) {
        try {
            return new Image(ImageDataFactory.create(path));
        } catch (MalformedURLException e) {
            throw new IOException("Could not read image from " + path);
        }
    }
}
