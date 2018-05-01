#include <dart/dart.hpp>
#include <dart/gui/gui.hpp>
#include <dart/utils/urdf/urdf.hpp>
#include <iostream>
#include <fstream>
#include <boost/circular_buffer.hpp>

using namespace std;
using namespace dart::common;
using namespace dart::dynamics;
using namespace dart::simulation;
using namespace dart::math;

class filter {
  public:
    filter(const int dim, const int n)
    {
      samples.set_capacity(n);
      total = Eigen::VectorXd::Zero(dim,1);
    }
    void AddSample(Eigen::VectorXd v)
    {
      if(samples.full()) 
      {
        total -= samples.front();
      }
      samples.push_back(v);
      total += v;
      average = total/samples.size();
    }
  
    boost::circular_buffer<Eigen::VectorXd> samples;
    Eigen::VectorXd total;
    Eigen::VectorXd average;
    
};

class MyWindow : public dart::gui::SimWindow
{
  public: 
    MyWindow(const WorldPtr& world)
    {
      setWorld(world);
      m3DOF = world->getSkeleton("m3DOF");
      qInit = m3DOF->getPositions();
      psi = 0; // Heading Angle
      steps = 0;
      outFile.open("constraints.csv");
      dqFilt = new filter(8, 100);
      cFilt = new filter(5, 100);
      R = 0.25;
      L = 0.68;//*6;
    }
    void timeStepping() override
    {
      // Read Positions, Speeds, Transform speeds to world coordinates and filter the speeds
      Eigen::Matrix<double, 4, 4> Tf = m3DOF->getBodyNode(0)->getTransform().matrix();
      psi =  atan2(Tf(0,0), -Tf(1,0));
      qBody1 = atan2(Tf(0,1)*cos(psi) + Tf(1,1)*sin(psi), Tf(2,1));
      Eigen::VectorXd q = m3DOF->getPositions();
      Eigen::VectorXd dq_orig = m3DOF->getVelocities();
      Eigen::Matrix<double, 8, 1> dq;
      dq << (Tf.block<3,3>(0,0) * dq_orig.head(3)) , (Tf.block<3,3>(0,0) * dq_orig.segment(3,3)), dq_orig(6), dq_orig(7);
      dqFilt->AddSample(dq);

      // Calculate the quantities we are interested in
      dpsi = dq(2);
      dpsiFilt = dqFilt->average(2);
      dqBody1 = -dq_orig(0);
      dqBody1Filt = (-dqFilt->average(0)*sin(psi) + dqFilt->average(1)*cos(psi));
      dthL = dq(6) + dqBody1;
      dthLFilt = dqFilt->average(6) + dqBody1Filt;
      dthR = dq(7) + dqBody1;
      dthRFilt = dqFilt->average(7) + dqBody1Filt;
      
      // Constraints
      // 1. dZ0 = 0                                               => dq_orig(4)*cos(qBody1) + dq_orig(5)*sin(qBody1) = 0
      // 2. da3 + R/L*(dthL - dthR) = 0                           => dq_orig(1)*cos(qBody1) + dq_orig(2)*sin(qBody1) + R/L*(dq_orig(6) - dq_orig(7)) = 0 
      // 3. da1*cos(psii) + da2*sin(psii) = 0                     => dq_orig(1)*sin(qBody1) - dq_orig(2)*cos(qBody1) = 0
      // 4. dX0*sin(psii) - dY0*cos(psii) = 0                     => dq_orig(3) = 0
      // 5. dX0*cos(psii) + dY0*sin(psii) - R/2*(dthL + dthR) = 0 => dq_orig(4)*sin(qBody1) - dq_orig(5)*cos(qBody1) - R/2*(dq_orig(6) + dq_orig(7) - 2*dq_orig(0)) = 0
      Eigen::Matrix<double, 5, 1> c;
      c << (dq_orig(4)*cos(qBody1) + dq_orig(5)*sin(qBody1)), (dq_orig(1)*cos(qBody1) + dq_orig(2)*sin(qBody1) + R/L*(dq_orig(6) - dq_orig(7))), (dq_orig(1)*sin(qBody1) - dq_orig(2)*cos(qBody1)), dq_orig(3), (dq_orig(4)*sin(qBody1) - dq_orig(5)*cos(qBody1) - R/2*(dq_orig(6) + dq_orig(7) - 2*dq_orig(0)));
      cFilt->AddSample(c);
      //if(Tf(2,1) > 0)
      {
      for(int i=0; i<8; i++) outFile << dq(i) << ", ";
      for(int i=0; i<8; i++) outFile << dqFilt->average(i) << ", ";
      outFile << psi << ", " << dpsi << ", " << qBody1 << ", " << dqBody1 << ", " <<  dthL << ", " << dthR << ", ";
      outFile << psiFilt << ", " << dpsiFilt << ", " << qBody1Filt << ", " << dqBody1Filt << ", " <<  dthLFilt << ", " << dthRFilt << ", ";
      for(int i=0; i<5; i++) outFile << c(i) << ", ";
      for(int i=0; i<5; i++) outFile << cFilt->average(i) << ", ";
      for(int i=0; i<8; i++) outFile << q(i) << ", "; 
      outFile << std::endl;
      }
      
      // arbitrary control inputs
      steps++;
      double headSign = (cos(2*3.14/80*steps*0.01) > 0) ? 1 : -1;
      double spinSign = (sin(2*3.14/20*steps*0.01) > 0) ? 1 : -1;
      double uL = ((steps < 1000) ? 0 : 5*headSign-5*spinSign);
      double uR = ((steps < 1000) ? 0 : 5*headSign+5*spinSign);
      /*double head = cos(2*3.14/10*steps*0.01);
      double spin = ( (sin(2*3.14/40*steps*0.01)>0) ? -1 : 1 );
      double uL = ((steps < 1000) ? 0 : 60*head-200*spin);
      double uR = ((steps < 1000) ? 0 : 60*head+200*spin);*/

      mForces << 0, 0, 0, 0, 0, 0, uL, uR;
      m3DOF->setForces(mForces);
      
      SimWindow::timeStepping();
    }
    ~MyWindow() {
      outFile.close();     
    }
    

  protected:

    SkeletonPtr m3DOF;

    Eigen::VectorXd qInit;

    Eigen::VectorXd dof1;

    double psi, dpsi, qBody1, dqBody1, dthL, dthR;
    double psiFilt, dpsiFilt, qBody1Filt, dqBody1Filt, dthLFilt, dthRFilt;

    double R;
    double L;
    
    int steps;

    Eigen::Matrix<double, 8, 1> mForces;
   
    ofstream outFile; 

    filter *dqFilt, *cFilt;
};


SkeletonPtr createFloor()
{
  SkeletonPtr floor = Skeleton::create("floor");

  // Give the floor a body
  BodyNodePtr body =
      floor->createJointAndBodyNodePair<WeldJoint>(nullptr).second;
//  body->setFrictionCoeff(1e16);

  // Give the body a shape
  double floor_width = 50;
  double floor_height = 0.05;
  std::shared_ptr<BoxShape> box(
        new BoxShape(Eigen::Vector3d(floor_width, floor_width, floor_height)));
  auto shapeNode
      = body->createShapeNodeWith<VisualAspect, CollisionAspect, DynamicsAspect>(box);
  shapeNode->getVisualAspect()->setColor(dart::Color::Blue());

  // Put the body into position
  Eigen::Isometry3d tf(Eigen::Isometry3d::Identity());
  tf.translation() = Eigen::Vector3d(0.0, 0.0, -floor_height / 2.0);
  body->getParentJoint()->setTransformFromParentBodyNode(tf);

  return floor;
}


SkeletonPtr create3DOF_URDF()
{
  // Load the Skeleton from a file
  dart::utils::DartLoader loader;
  SkeletonPtr threeDOF = 
      loader.parseSkeleton("/home/krang/dart/09-URDF/3DOF-WIP/3dof.urdf");
      //loader.parseSkeleton("/home/panda/myfolder/wholebodycontrol/09-URDF/3DOF-WIP/3dof.urdf");
  threeDOF->setName("m3DOF");

  // Set the body inertia to a sensible value
  SkeletonPtr krangFixedWheel =
      loader.parseSkeleton("/home/krang/dart/09-URDF/KrangFixedWheels/krang_fixed_wheel.urdf");
  krangFixedWheel->setName("m18DOF");
  Eigen::Matrix<double, 18, 1> qInit;
  qInit << -M_PI/4, -4.588, 0.0, 0.0, 0.0548, -1.0253, 0.0, -2.1244, -1.0472, 1.5671, 0.0, -0.0548, 1.0253, 0.0, 2.1244, 1.0472, 0.0037, 0.0;
  krangFixedWheel->setPositions(qInit);
  int nBodies = krangFixedWheel->getNumBodyNodes();
  for(int i=0; i<nBodies; i++){
    dart::dynamics::BodyNode* body= krangFixedWheel->getBodyNode(i);
    std::cout << body->getName() << std::endl;
  }

  
  // Get it into a useful configuration
  double psiInit = M_PI/4, qBody1Init = M_PI/6;
  Eigen::Transform<double, 3, Eigen::Affine> baseTf = Eigen::Transform<double, 3, Eigen::Affine>::Identity();
  // RotX(pi/2)*RotY(-pi/2+psi)*RotX(-qBody1)
  baseTf.prerotate(Eigen::AngleAxisd(-qBody1Init,Eigen::Vector3d::UnitX())).prerotate(Eigen::AngleAxisd(-M_PI/2+psiInit,Eigen::Vector3d::UnitY())).prerotate(Eigen::AngleAxisd(M_PI/2, Eigen::Vector3d::UnitX()));
  Eigen::AngleAxisd aa(baseTf.matrix().block<3,3>(0,0));
  Eigen::Matrix<double, 8, 1> q;
//  q << 1.2092, -1.2092, -1.2092, 0, 0, 0.28, 0, 0;
  q << aa.angle()*aa.axis(), 0, 0, 0.28, 0, 0;
  threeDOF->setPositions(q);

  return threeDOF;
}


int main(int argc, char* argv[])
{

  SkeletonPtr threeDOF = create3DOF_URDF();
  SkeletonPtr floor = createFloor();

  WorldPtr world = std::make_shared<World>();
  world->addSkeleton(threeDOF);
  world->addSkeleton(floor);

  MyWindow window(world);
  glutInit(&argc, argv);
  window.initWindow(1280,720, "3DOF URDF");
  glutMainLoop();
}
