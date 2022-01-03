
plugins {
    id("com.hartwig.java-conventions")
    application
}

description = "HMF Tools - Sigs"

dependencies {
    implementation(project(":hmf-common"))
    implementation(project(":patient-db"))
    testImplementation(libs.junit)
}

application {
    mainClass.set("com.hartwig.hmftools.sigs.fitter.SampleFitter")
}
