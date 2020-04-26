if [ $# -eq 1 ]
then
    echo -e "\e[31mDeleting build directory\e[0m"
    rm -r build
    echo -e "\e[31mMaking build directory\e[0m"
    mkdir build && cd build
    echo -e "\e[31mCalling CMake\e[0m"
    cmake ..
else
    cd build
fi
echo -e "\e[31mCalling make\e[0m"
make
echo -e "\e[31mRunning project...\e[0m"
./project -ksp_type gmres -pc_type lu
echo -e "\e[31mDone!\e[0m"
