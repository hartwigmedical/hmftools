package com.hartwig.hmftools.patientreporter.cfreport.components;

import java.net.MalformedURLException;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.itextpdf.io.IOException;
import com.itextpdf.io.image.ImageDataFactory;
import com.itextpdf.layout.element.Div;
import com.itextpdf.layout.element.Image;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Text;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ReportSignature {

    private ReportSignature() {
    }

    @NotNull
    public static Paragraph createEndOfReportIndication() {
        return new Paragraph("— End of report —").setMarginTop(40).addStyle(ReportResources.smallBodyTextStyle());
    }

    @NotNull
    public static Div createSignatureDiv(@NotNull String rvaLogoPath, @NotNull String signaturePath) throws IOException {
        Div div = new Div();
        div.setKeepTogether(true);
        div.setMarginTop(40);

        try {
            Image rvaLogo = new Image(ImageDataFactory.create(rvaLogoPath));
            rvaLogo.setMaxHeight(58);
            div.add(rvaLogo);
        } catch (MalformedURLException e) {
            throw new IOException("Failed to read RVA logo image at " + rvaLogoPath);
        }

        Paragraph signatureText =
                new Paragraph().setFont(ReportResources.fontBold()).setFontSize(10).setFontColor(ReportResources.PALETTE_BLACK);

        signatureText.add(ReportResources.SIGNATURE_NAME + ",\n");
        signatureText.add(new Text(ReportResources.SIGNATURE_TITLE).setFont(ReportResources.fontRegular()));
        div.add(signatureText);

        try {
            Image signatureImage = new Image(ImageDataFactory.create(signaturePath));
            signatureImage.setMaxHeight(60);
            signatureImage.setMarginTop(-20); // Set negative margin so the signature slightly overlaps the signature text
            signatureImage.setMarginLeft(10);
            div.add(signatureImage);
        } catch (MalformedURLException e) {
            throw new IOException("Failed to read signature image at " + signaturePath);
        }

        return div;
    }

    @NotNull
    public static Div createSignatureDivPanel(@NotNull String signaturePath) throws IOException {
        Div div = new Div();
        div.setKeepTogether(true);
        div.setMarginTop(40);

        Paragraph signatureText =
                new Paragraph().setFont(ReportResources.fontBold()).setFontSize(10).setFontColor(ReportResources.PALETTE_BLACK);

        signatureText.add(ReportResources.SIGNATURE_NAME + ",\n");
        signatureText.add(new Text(ReportResources.SIGNATURE_TITLE).setFont(ReportResources.fontRegular()));
        div.add(signatureText);

        try {
            Image signatureImage = new Image(ImageDataFactory.create(signaturePath));
            signatureImage.setMaxHeight(60);
            signatureImage.setMarginTop(-20); // Set negative margin so the signature slightly overlaps the signature text
            signatureImage.setMarginLeft(10);
            div.add(signatureImage);
        } catch (MalformedURLException e) {
            throw new IOException("Failed to read signature image at " + signaturePath);
        }

        return div;
    }
}
