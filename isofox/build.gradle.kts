
plugins {
    id("com.hartwig.java-conventions")
    application
}

description = "HMF Tools - Isofox"

dependencies {
    implementation(project(":hmf-common"))
    implementation(project(":patient-db"))
    testImplementation(libs.junit)
}

application {
    mainClass.set("com.hartwig.hmftools.isofox.Isofox")
}
