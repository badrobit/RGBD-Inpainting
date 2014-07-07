echo "Compiling & Installing RGB-D Inpainting"
mkdir build
cd build
scriptpath="PWD"
echo "Installing any required compilation tools" 
sudo apt-get -qq install cmake cmake-curses-gui g++
opencvthere=$(pkg-config --libs opencv | grep \".*libopencv.*\")
if [ -z $opencvthere ]
then
	echo "OpenCV Found Skipping Install"
else
	cd ~
	if [ ! -e OpenCV.sh ]
	then
		wget https://gist.githubusercontent.com/Bad-Robot/f4a28666fc14f1af06d5/raw/f6ca8605b7874b2ac9bd6bbd8e5a2500e9ea04c1/OpenCV.sh
	fi
	chmod a+x OpenCV.sh
	sudo sh OpenCV.sh
	cd $scriptpath
fi

echo "Running CMake ..."
cmake ..
make -j4
clear
echo "RGB-D Inpainting Compiled and Installed"
echo "From now on you can use cmake and make commands alone."
echo "If you make changes to the source please push them back to the original repository!"
