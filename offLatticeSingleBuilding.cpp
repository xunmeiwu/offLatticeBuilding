#include "palabos3D.h"
#include "palabos3D.hh"
#include <vector>
#include <cmath>
#include <iostream>
//#include <sstream>
//#include <fstream>
#include <string>
#include <iomanip>

using namespace plb;
using namespace plb::descriptors;
using namespace std;

typedef double T;
typedef Array<T, 3> Velocity;
#define DESCRIPTOR descriptors::MRTD3Q19Descriptor

std::string outputDir("./tmp/");

static T SMALL = 1.0e-30;//doubleScalar.H
T aveFrequncy = SMALL; //This keeps the sampling amount for the time average

/** 案例参数 **/
std::string geometry_fname = "triangle.stl";  //读取的模型名字
std::string caseNameChar = "offLatticeSingleBuilding"; //案例名称，log文件及vtk文件保存名.

static plint xDirection = 0;
static plint borderWidth = 1;  // Because Guo acts in a one-cell layer.
// Requirement: margin>=borderWidth.
static plint margin =
        1;  // Extra margin of allocated cells around the obstacle, for the case of moving walls.
static plint blockSize = 0;  // Size of blocks in the sparse/parallel representation.
// Zero means: don't use sparse representation.
static plint extendedEnvelopeWidth = 2;  // Extrapolated off-lattice BCs.

const T X_Length = 4.08;            //流域X方向长度；
const T Y_Length = 2.2;            //流域Y方向长度；
const T Z_Length = 1.8;            //流域Z方向长度；
const T cx = 1.62;                 //Position of the obstacle, in physical units.
const T cy = 1.1;
const T stlLength = 0.16;         //stl模型的长度

/** 单位转换参数 **/
const T resolution = 25.;                   //物理长度的网格数
const T charNu = 1.5e-5;                //代表物理动粘性系数；
const T latticeU = 0.01;        //晶格速度，与马赫数匹配

const T rho0 = 1.0;
const Array<T, 3> u0((T) 0., (T) 0., (T) 0.);


/** 模拟参数 **/
const T iniT = (T) 0;                 // 模拟初始时刻，单位s，一般从0开始，若大于0则需要读取现有晶格
const T maxT = (T) iniT + 3. + 0.1; // 模拟截止时刻，单位s
const T cSmago = 0.12;              // Smagorisky常数

const T statT = (T) 0.1;     // 控制台Log显示时间间隔，单位s
const T vtkStartT = (T) 0.; //VTK文件输出开始时间，单位s
const T vtkT = (T) 1;    // VTK文件输出时间间隔，单位s
const T imSave = (T) 0.1;    // image文件输出时间间隔，单位s
const bool useAve = true;    //是否进行时间平均
const T aveStartT = (T) 2.; //平均开始时间
const T aveDelta = (T) 0.05; //平均时间间隔，单位s

const std::string inletPath = "samplingFace/"; // Path and prefix of inlet files
const T inletFileScaling = 1.;
const Array<T, 3> ratioPntToLattice(1., 1.2 / 2.2, 1.0 / 1.8); //Size ratio of two point set, =FOAM/lattice
const Array<T, 3> readPntStartPosition(-0.5, -0.6, 0.);        //Start position of FOAM points for lattice
const Array<plint, 3> pntDim(1, 122, 142);                     //Dimension of the inlet FOAM points

const bool latSave = false;                                                  //是否保存lattice文件
const T latSaveT = (T) 60;                                                  // lattice保存时间间隔，单位s
const bool latLoad = false;                                                  //是否读取lattice文件
const T latLoadT = 1000;                                                    // lattice读取开始时间，单位s
std::string latSaveFile = "checkpointing_", latLoadFile = "checkpointing_"; //写入与读取lattice文件的文件名


/** Get formated current system time and **/
string getTime() {
    time_t timep;
    time(&timep);
    char tmp[64];
    strftime(tmp, sizeof(tmp), "%Y-%m-%d %H:%M:%S", localtime(&timep));
    return tmp;
}

//#include "readInlet_FOAM.h"

/** 确定模型每个边的边界条件 **/
void boundarySetup(MultiBlockLattice3D<T, DESCRIPTOR> &lattice,
                   IncomprFlowParam<T> parameters,
                   OnLatticeBoundaryCondition3D<T, DESCRIPTOR> &boundaryCondition) {
    Array<T, 3> uBoundary(5 * latticeU, (T) 0., (T) 0.);
    plint nx = parameters.getNx();
    plint ny = parameters.getNy();
    plint nz = parameters.getNz();

    Box3D Roof(0, nx - 1, 0, ny - 1, nz - 1, nz - 1);        //Roof
    Box3D Ground(0, nx - 1, 0, ny - 1, 0, 0);              //Ground
    Box3D Inlet(0, 0, 0, ny - 1, 0, nz - 1);                //Leftside 入口处的四边交角归纳到Roof，Ground等
    Box3D Outlet(nx - 1, nx - 1, 0, ny - 1, 0, nz - 1);    //Rightside 出口处的四边交角归纳到Roof，Ground等
    Box3D Leftside(0, nx - 1, 0, 0, 0, nz - 1);            //Frontside
    Box3D Rightside(0, nx - 1, ny - 1, ny - 1, 0, nz - 1);    //Backside

    //入口dirichlet速度条件
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, Roof, boundary::freeslip);
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, Ground, boundary::dirichlet);
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, Inlet, boundary::dirichlet);
    //出口pressure条件，0P为x轴正方向，1N为Y轴负方向
    boundaryCondition.addPressureBoundary0P(Outlet, lattice, boundary::dirichlet);
//    出口outflow条件(zero velocity-gradient)
//    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, Outlet, boundary::outflow);
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, Leftside, boundary::dirichlet);
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, Rightside, boundary::dirichlet);

    pcout << "Roof:(" << Roof.x0 << "," << Roof.x1 << "," << Roof.y0 << "," << Roof.y1 << "," << Roof.z0 << ","
          << Roof.z1 << ")." << endl;
    pcout << "Ground:(" << Ground.x0 << "," << Ground.x1 << "," << Ground.y0 << "," << Ground.y1 << "," << Ground.z0
          << "," << Ground.z1 << ")." << endl;
    pcout << "Inlet:(" << Inlet.x0 << "," << Inlet.x1 << "," << Inlet.y0 << "," << Inlet.y1 << "," << Inlet.z0 << ","
          << Inlet.z1 << ")." << endl;
    pcout << "Outlet:(" << Outlet.x0 << "," << Outlet.x1 << "," << Outlet.y0 << "," << Outlet.y1 << "," << Outlet.z0
          << "," << Outlet.z1 << ")." << endl;
    pcout << "Leftside:(" << Leftside.x0 << "," << Leftside.x1 << "," << Leftside.y0 << "," << Leftside.y1 << ","
          << Leftside.z0 << "," << Leftside.z1 << ")." << endl;
    pcout << "Rightside:(" << Rightside.x0 << "," << Rightside.x1 << "," << Rightside.y0 << "," << Rightside.y1 << ","
          << Rightside.z0 << "," << Rightside.z1 << ")." << endl;

    //所有位置给定初始速度和密度并初始化
    setBoundaryVelocity(lattice, lattice.getBoundingBox(), u0);
    setBoundaryDensity(lattice, Outlet, rho0);
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), rho0, u0);
    setBoundaryVelocity(lattice, Inlet, uBoundary);
    lattice.initialize();
}

/** VTK file output **/
void writeVTK(MultiBlockLattice3D<T, DESCRIPTOR> &lattice,
              IncomprFlowParam<T> const &parameters,
              plint iter,
              MultiTensorField3D<T, 3> &velocity,
              MultiScalarField3D<T> &rho,
              MultiScalarField3D<T> &nut,
              MultiTensorField3D<T, 3> &uSum,
              MultiTensorField3D<T, 6> &uPrimeSum,
              MultiScalarField3D<T> &nutSum,
              MultiScalarField3D<T> &nutSquareSum,
              MultiScalarField3D<T> &rhoSum,
              MultiScalarField3D<T> &rhoSquareSum,
              MultiTensorField3D<T, 3> &uSquareSum
) {
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();
    T pressureScale = (dx * dx) / (dt * dt) * DESCRIPTOR<T>::cs2;

    //输出各种瞬时值到vtk文件
    //Obtain the iteration time
    std::ostringstream strOS;
    strOS << fixed << setprecision(6) << iter * dt;
    std::string strIterTime = strOS.str();
    VtkImageOutput3D<T> vtkOut("vtk" + strIterTime + "s", dx);
    vtkOut.writeData<3, T>(velocity, "velocity", dx / dt);
    vtkOut.writeData<T>(rho, "rho", 1.);
    vtkOut.writeData<T>(rho, "pressure", pressureScale);
    vtkOut.writeData<T>(nut, "nut", util::sqr(dx) / dt);

    //Calculate meanVelocity
    std::unique_ptr<MultiTensorField3D<T, 3>> uMean = multiply((T) 1. / aveFrequncy, uSum,
                                                               lattice.getBoundingBox());

    //Calculate uPrimeInstantaneous <u'u'><v'v'><w'w'><u'v'><u'w'><v'w'>
    std::unique_ptr<MultiTensorField3D<T, 6>> uPrimeInstantaneous = computeInstantaneousReynoldsStress(velocity, *uMean,
                                                                                                       lattice.getBoundingBox());
    vtkOut.writeData<6, T>(*uPrimeInstantaneous, "uPrimeInstantaneous", util::sqr(dx / dt));

    //**************************************************************************************

    if (useAve && iter >= parameters.nStep(aveStartT)) {
        //Calculate uPrime2Mean <u'u'><v'v'><w'w'><u'v'><u'w'><v'w'>
        std::unique_ptr<MultiTensorField3D<T, 6>> uPrime2Mean = multiply((T) 1. / aveFrequncy, uPrimeSum,
                                                                         lattice.getBoundingBox());

        //Calculate meanNut <nut>
        std::unique_ptr<MultiScalarField3D<T>> nutMean = multiply((T) 1. / aveFrequncy, nutSum,
                                                                  lattice.getBoundingBox());

        //Calculate nutPrime2Mean, <Nut'^2>=<Nut^2>-<Nut>^2
        std::unique_ptr<MultiScalarField3D<T>> nutPrime2Mean = subtract(
                *multiply((T) 1. / aveFrequncy, nutSquareSum, lattice.getBoundingBox()),
                *multiply(*nutMean, *nutMean, lattice.getBoundingBox()),
                lattice.getBoundingBox());

        //Calculate rhoMean <rho>
        std::unique_ptr<MultiScalarField3D<T>> rhoMean = multiply((T) 1. / aveFrequncy, rhoSum,
                                                                  lattice.getBoundingBox());

        //Calculate rhoPrime2Mean <Rho'^2>=<Rho^2>-<Rho>^2
        std::unique_ptr<MultiScalarField3D<T>> rhoPrime2Mean = subtract(
                *multiply((T) 1. / aveFrequncy, rhoSquareSum, lattice.getBoundingBox()),

                *multiply(*rhoMean, *rhoMean, lattice.getBoundingBox()),
                lattice.getBoundingBox());
        //Calculate ((X1^2+X2^2+X3^2+...+Xn^2)/n)^(1/2)
        std::unique_ptr<MultiTensorField3D<T, 3>> RMS = computeSqrt(
                *multiply((T) 1. / aveFrequncy, uSquareSum, lattice.getBoundingBox()),
                lattice.getBoundingBox());

        //VTK output
        vtkOut.writeData<3, T>(*uMean, "uMean", dx / dt);
        vtkOut.writeData<6, T>(*uPrime2Mean, "uPrime2Mean", util::sqr(dx / dt));
        vtkOut.writeData<T>(*rhoMean, "pMean", pressureScale);
        vtkOut.writeData<T>(*rhoPrime2Mean, "pPrime2Mean", pressureScale * util::sqr(dx / dt));
        vtkOut.writeData<T>(*nutMean, "nutMean", util::sqr(dx) / dt);
        vtkOut.writeData<T>(*nutPrime2Mean, "nutPrime2Mean", util::sqr(util::sqr(dx) / dt));
        vtkOut.writeData<T>(*rhoMean, "rhoMean", 1.);
        vtkOut.writeData<T>(*rhoPrime2Mean, "rhoPrime2Mean", 1.);
        vtkOut.writeData<3, T>(*RMS, "Urms", dx / dt);

    }
}

void runProgram(IncomprFlowParam<T> parameters) {
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();
    plint nx = parameters.getNx();
    plint ny = parameters.getNy();
    plint nz = parameters.getNz();

    Box3D boxLB0(0, nx - 1, 0, ny - 1, 0, nz - 1);

    /*
     * 读取stl文件
     */
    pcout << std::endl << "Reading STL data for the obstacle geometry." << std::endl;
    // The triangle-set defines the surface of the geometry.
    TriangleSet<T> *triangleSet = new TriangleSet<T>(geometry_fname, DBL);
//    triangleSet->scale(parameters.getDeltaX()); // In lattice units from now on...
    Cuboid<T> bCuboid = triangleSet->getBoundingCuboid();
    Array<T, 3> obstacleCenter = (T) 0.5 * (bCuboid.lowerLeftCorner + bCuboid.upperRightCorner);
    triangleSet->translate(-obstacleCenter);
    triangleSet->scale(stlLength * resolution);
    Array<T, 3> cornerLB(cx * resolution,
                         cy * resolution,
                         stlLength * resolution + 1);
    triangleSet->translate(cornerLB);
    triangleSet->writeBinarySTL(outputDir + "obstacle_LB.stl");


    pcout << "stlLength * resolution =" << stlLength * resolution << "  cx * resolution =" << cx * resolution
          << "  cy * resolution =" << cy * resolution << std::endl;

    DEFscaledMesh<T> *defMesh = new DEFscaledMesh<T>(*triangleSet, 0, xDirection,
                                                     margin,
                                                     Dot3D(0, 0, 0));

    delete triangleSet;
    triangleSet = 0;
    TriangleBoundary3D<T> triangleBoundary(*defMesh);
    delete defMesh;
    defMesh = 0;
    triangleBoundary.getMesh().inflate();

    /*
     * 像素化计算区域
     */
    pcout << std::endl << "Voxelizing the domain." << std::endl;
    const int flowType = voxelFlag::outside;
    VoxelizedDomain3D<T> voxelizedDomain(
            triangleBoundary, flowType, boxLB0, borderWidth, extendedEnvelopeWidth,
            blockSize);
    pcout << getMultiBlockInfo(voxelizedDomain.getVoxelMatrix()) << std::endl;
    {
        VtkImageOutput3D<T> vtkOut(outputDir + "voxels_full_domain", 1);
        vtkOut.writeData<float>(*copyConvert<int, T>(voxelizedDomain.getVoxelMatrix(), boxLB0), "voxel", 1.0);
    }

    /*
     * Generate the lattice, the density and momentum blocks.
     */

    //定义lattice,相当于MultiBlockLattice3D<T, DESCRIPTOR> *lattice = new MultiBlockLattice3D<T, DESCRIPTOR>
    MultiBlockLattice3D<T, DESCRIPTOR> *lattice =
            new MultiBlockLattice3D<T, DESCRIPTOR>(voxelizedDomain.getVoxelMatrix());
    //定义外部流场边界条件
    defineDynamics(
            *lattice, lattice->getBoundingBox(),
            new SmagorinskyMRTdynamics<T, DESCRIPTOR>(parameters.getOmega(), cSmago));
    pcout << "Using MRT dynamics." << std::endl;
    //定义内部流场边界条件noDynamic
    defineDynamics(
            *lattice, voxelizedDomain.getVoxelMatrix(), lattice->getBoundingBox(),
            new NoDynamics<T, DESCRIPTOR>(), voxelFlag::inside);
    defineDynamics(
            *lattice, voxelizedDomain.getVoxelMatrix(), lattice->getBoundingBox(),
            new NoDynamics<T, DESCRIPTOR>(), voxelFlag::innerBorder);
    defineDynamics(
            *lattice, voxelizedDomain.getVoxelMatrix(), lattice->getBoundingBox(),
            new NoDynamics<T, DESCRIPTOR>(), voxelFlag::undetermined);
    defineDynamics(
            *lattice, voxelizedDomain.getVoxelMatrix(), lattice->getBoundingBox(),
            new NoDynamics<T, DESCRIPTOR>(), voxelFlag::toBeInside);


    /*
     * 生成stl的off-lattice边界条件和outer边界条件。
     */

    pcout << "Generating boundary conditions." << std::endl;

    OffLatticeBoundaryCondition3D<T, DESCRIPTOR, Velocity> *boundaryCondition;

    BoundaryProfiles3D<T, Velocity> profiles;
    bool useAllDirections = true;//Extrapolation scheme for the off lattice boundary condition.
    OffLatticeModel3D<T, Velocity> *offLatticeModel = 0;
    profiles.setWallProfile(new NoSlipProfile3D<T>);
    offLatticeModel = new GuoOffLatticeModel3D<T, DESCRIPTOR>(
            new TriangleFlowShape3D<T, Array<T, 3> >(voxelizedDomain.getBoundary(), profiles), flowType,
            useAllDirections);
//    offLatticeModel = new FilippovaHaenelLocalModel3D<T,DESCRIPTOR>(new TriangleFlowShape3D<T, Array<T, 3> >(voxelizedDomain.getBoundary(), profiles), flowType);
//    bool velIsJ = true;
//    offLatticeModel->setVelIsJ(velIsJ);
    boundaryCondition = new OffLatticeBoundaryCondition3D<T, DESCRIPTOR, Velocity>(
            offLatticeModel, voxelizedDomain, *lattice);

    boundaryCondition->insert();

    //The boundary condition algorithm or the outer domain.
    OnLatticeBoundaryCondition3D<T, DESCRIPTOR> *outerBoundaryCondition =
            createInterpBoundaryCondition3D<T, DESCRIPTOR>();
    pcout << "InterpolationBC3D...OK" << endl;
    boundarySetup(*lattice, parameters, *outerBoundaryCondition);
    delete outerBoundaryCondition;


    writeLogFile(parameters, caseNameChar); //输出转换参数
    pcout << "Omega= " << parameters.getOmega() << std::endl;
    pcout << "Tao= " << parameters.getTau() << std::endl;
    pcout << "latticeU= " << latticeU << std::endl;
    pcout << "getRe= " << parameters.getRe() << std::endl;
    pcout << "Resolution= " << parameters.getResolution() << std::endl;
    pcout << "DeltaX= " << parameters.getDeltaX() << std::endl;
    pcout << "DeltaT= " << parameters.getDeltaT() << std::endl;
    pcout << "Nx=" << parameters.getNx() << ", Ny=" << parameters.getNy() << ", Nz=" << parameters.getNz() << std::endl;

//    MultiBlockLattice3D<T, DESCRIPTOR> *lattice = new MultiBlockLattice3D<T, DESCRIPTOR>(parameters.getNx(),
//                                                                                         parameters.getNy(),
//                                                                                         parameters.getNz(),
//                                                                                         new SmagorinskyMRTdynamics<T, DESCRIPTOR>(
//                                                                                                 parameters.getOmega(),
//                                                                                                 cSmago));

    //define basic instantaneous physical quantities
    MultiTensorField3D<T, 3> velocity(nx, ny, nz, u0); //velocity
    MultiScalarField3D<T> nut(nx, ny, nz, (T) 0.);                 //nut
    MultiScalarField3D<T> rho(nx, ny, nz, (T) 0.);                 //rho
    //Define all tensor and scalar fields for averaged operation
    MultiTensorField3D<T, 3> uSum(nx, ny, nz, u0);
    MultiTensorField3D<T, 6> uPrimeSum(nx, ny, nz);
    MultiScalarField3D<T> nutSum(nx, ny, nz, (T) 0.);
    MultiScalarField3D<T> nutSquareSum(nx, ny, nz, (T) 0.);
    MultiScalarField3D<T> rhoSum(nx, ny, nz, (T) 0.);
    MultiScalarField3D<T> rhoSquareSum(nx, ny, nz, (T) 0.);
    MultiTensorField3D<T, 3> uSquareSum(nx, ny, nz, u0);

////定义边界条件
//    OnLatticeBoundaryCondition3D<T, DESCRIPTOR> *boundaryCondition;
//    pcout << "Setting BoundaryCondition:";
//    boundaryCondition = createInterpBoundaryCondition3D<T, DESCRIPTOR>();
//    pcout << "InterpolationBC3D...OK" << endl;
//    boundarySetup(*lattice, parameters,*boundaryCondition);
//    delete boundaryCondition;

    std::ostringstream strOS;

    // 读取之前的晶格数据
    if (latLoad) {
        //loadRawMultiBlock(lattice, "checkpoint.dat");
        strOS << latLoadT;
        std::string strLoadTime = strOS.str();
        pcout << "Loading lattice file " << latLoadFile + strLoadTime + "s.dat"
              << " ...";
        loadBinaryBlock(*lattice, latLoadFile + strLoadTime + "s.dat");
        pcout << "OK." << endl;
    }

    global::timer("mainLoop").start(); //主循环计时器
    global::timer("iteLog").start();   // log显示之间的循环计时器

//    FOAM_InletBC InletBCFromFOAM(inletPath, Box3D(0, 0, 1, ny - 2, 1, nz - 2), inletFileScaling, pntDim,
//                                 ratioPntToLattice, readPntStartPosition);
//    InletBCFromFOAM.initializeFOAMPnts(parameters);

    // 从初始时间iniT开始进行主循环
    pcout << "Starting iteration, Current time: " << getTime() << endl;
    for (plint iT = iniT / dt; iT * dt < maxT; ++iT) {

        //Update new basic physical quantities
        if (iT % parameters.nStep(statT) == 0 ||
            (useAve && iT % parameters.nStep(aveDelta) == 0 && iT >= parameters.nStep(aveStartT) &&
             iT > parameters.nStep(iniT)) ||
            (iT % parameters.nStep(vtkT) == 0 && iT >= parameters.nStep(vtkStartT) && iT > parameters.nStep(iniT)) ||
            iT % parameters.nStep(imSave) == 0) {
            rho = *computeDensity(*lattice, lattice->getBoundingBox());
            velocity = *computeVelocity(*lattice, lattice->getBoundingBox());
        }

        /// 输出控制台信息
        if (iT % parameters.nStep(statT) == 0) {
            pcout << endl;
            pcout << "step " << iT << "; t=" << iT * dt << endl;
            Array<T, 3> nearGroundVel;
            Array<T, 3> nearRoofVel;
            for (plint i = 0; i < 3; ++i) {
                nearGroundVel[i] = (*extractComponent(velocity, i)).get(nx / 2, ny / 2, 1);
                nearRoofVel[i] = (*extractComponent(velocity, i)).get(nx / 2, ny / 2, nz - 2);
            }

            pcout << "nearGroundVel= " << fixed << setprecision(6) << nearGroundVel[0] * dx / dt << ","
                  << nearGroundVel[1] * dx / dt << "," << nearGroundVel[2] * dx / dt << ";";
            pcout << "nearRoofVel= " << fixed << setprecision(6) << nearRoofVel[0] * dx / dt << ","
                  << nearRoofVel[1] * dx / dt << "," << nearRoofVel[2] * dx / dt << std::endl;
        }

//        //Reading FOAM inlet U data
//        InletBCFromFOAM.readFOAMU(*lattice, parameters, iT);

        lattice->collideAndStream();

        ///输出image
        if (iT % parameters.nStep(imSave) == 0) {
            pcout << "Writing Gif ...";
            const plint imSize = 600;
            Box3D slice(0, nx - 1, ny / 2, ny / 2, 0, nz - 1);
            ImageWriter<T> imageWriter("leeloo");
            imageWriter.writeScaledGif(createFileName("ux", iT, 6),
                                       *computeVelocityComponent(*lattice, slice, 0),
                                       imSize, imSize);
            imageWriter.writeScaledGif(createFileName("uz", iT, 6),
                                       *computeVelocityComponent(*lattice, slice, 2),
                                       imSize, imSize);
            imageWriter.writeScaledGif(createFileName("velNorm", iT, 6),
                                       *computeVelocityNorm(*lattice, slice),
                                       imSize, imSize);
            pcout << "OK." << endl;
        }

        ///输出VTK文件
        if (iT % parameters.nStep(vtkT) == 0 && iT >= parameters.nStep(vtkStartT) && iT > parameters.nStep(iniT)) {
            pcout << "Saving VTK file ...";
            writeVTK(*lattice, parameters, iT, velocity, rho, nut, uSum, uPrimeSum, nutSum, nutSquareSum, rhoSum,
                     rhoSquareSum,
                     uSquareSum);
            pcout << "OK." << endl;
        }

        ///计算平均值
        if (useAve && iT % parameters.nStep(std::max(dt, aveDelta)) == 0 &&
            iT >= parameters.nStep(aveStartT) &&
            iT > parameters.nStep(iniT)) //iT>parameters.nStep(iniT)是必须的，防止初始读入数据0参与平均
        {
            if (iT % parameters.nStep(statT) == 0) {
                pcout << "Calculate Average Velocity ...";
            }
            aveFrequncy = aveFrequncy + 1;
            //Sum instantaneous velocity into uSum
            add(uSum, velocity, uSum, lattice->getBoundingBox());
            //Calculate uMean
            std::unique_ptr<MultiTensorField3D<T, 3>> uMean = multiply((T) 1. / aveFrequncy, uSum,
                                                                       lattice->getBoundingBox());
            //Calculate instantaneous UPrime u'u', v'v', w'w', u'v', u'w', v'w'
            std::unique_ptr<MultiTensorField3D<T, 6>> currentUPrime = computeInstantaneousReynoldsStress(velocity,
                                                                                                         *uMean,
                                                                                                         lattice->getBoundingBox());
            //Sum instantaneous UPrime into uPrimeSum
            addInPlace(uPrimeSum, *currentUPrime);
            //Sum instantaneous nut into nutSum
            addInPlace(nutSum, nut);
            //Sum nutSquare into nutSquareSum
            addInPlace(nutSquareSum, *multiply(nut, nut, lattice->getBoundingBox()));
            //Sum rho into rhoSum
            addInPlace(rhoSum, rho);
            //Sum rhoSquare into rhoSquareSum
            addInPlace(rhoSquareSum, *multiply(rho, rho, lattice->getBoundingBox()));
            if (iT % parameters.nStep(statT) == 0) {
                pcout << "OK" << endl;
            }
        }

        ///保存晶格状态
        if (latSave && iT % parameters.nStep(latSaveT) == 0 && iT > parameters.nStep(iniT)) {
            //saveRawMultiBlock(lattice, "checkpoint.dat");
            std::ostringstream strOS1;
            strOS1 << fixed << setprecision(6) << iT * dt;
            //strIterTime = strOS1.str();
            std::string strLatticeSave = latSaveFile + strOS1.str() + "s.dat";
            pcout << "Saving the lattice file: " << strLatticeSave << " ...";
            saveBinaryBlock(*lattice, strLatticeSave);
            pcout << "OK." << endl;
        }

        if (iT % parameters.nStep(statT) == 0) {
            pcout << "av energy="
                  << setprecision(10) << getStoredAverageEnergy<T>(*lattice)
                  << "; av rho="
                  << setprecision(10) << getStoredAverageDensity<T>(*lattice) << endl;
            pcout << "Time spent during previous iterations: "
                  << global::timer("iteLog").stop()
                  << ". Current time: " << getTime() << endl;
            global::timer("iteLog").restart(); // 重启log显示之间的循环计时器
        }
    }
    pcout << "Iteration finished." << endl;
    pcout << "Total time spent during iterations: "
          << global::timer("mainLoop").stop() << endl;
    pcout << "Current time: " << getTime();

    delete lattice;
}


/** 计算主程序 **/
int main(int argc, char *argv[]) {
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");
//确定转换参数
    IncomprFlowParam<T> parameters(
            (T) latticeU,            // Re
            (T) 66666.7,
            resolution, // N,Resolution
            X_Length, // lx
            Y_Length, // ly
            Z_Length  // lz
    );

    runProgram(parameters);
}
